TcGSA.LR.alt <-
function(expr, gmt, Patient_ID, TimePoint, func = "linear", minGSsize=10, maxGSsize=500, group.var=NULL, separatePatients=FALSE){
  
  #DEBUG: data(data_simu_TcGSA);
  #expr=expr; Patient_ID=design$Patient_ID; TimePoint=design$Timepoint; gmt=gmt_sim; func="linear"; gs=2; maxGSsize=500; group.var=NULL;
  #for (i in which(tcGSA_cubic_modules$CVG_H1>6)){
  #  print(paste(i, ":", length(intersect(gmt$genesets[[i]], rownames(expr)))))
  #}
  
  if(!is.null(group.var) && separatePatients){
    stop("'separatePatients' is TRUE while 'group.var' is not NULL.\n This is an attempt to separate patients in a multiple group setting.\n This is not handled by the TcGSA.LR function.\n\n")
  }

#   library(GSA)
#   library(lme4)
#   library(reshape2)
#   require(splines)
   
  LR <- numeric(length(gmt$genesets))
  CVG_H0 <- numeric(length(gmt$genesets))
  CVG_H1 <- numeric(length(gmt$genesets))
  estim_expr <- list()
  splines_DF <- rep(NA, length(gmt$genesets))
  
  for (gs in 1:length(gmt$genesets)){
    probes <- intersect(gmt$genesets[[gs]], rownames(expr))

    
    if(length(probes)>0 && length(probes)<=maxGSsize && length(probes)>=minGSsize){                                                       
      expr_temp <- t(expr[probes, ])
      rownames(expr_temp) <- NULL
      
      if(is.null(group.var)){
        data_temp <- cbind.data.frame(Patient_ID, TimePoint, expr_temp)
        data_lm <- melt(data_temp, id.vars=c("Patient_ID", "TimePoint"), variable.name ="probe", value.name="expression")
        colnames(data_lm)[2] <- "t1"
        
        if(func=="cubic"){
          data_lm$t2 <- (data_lm$t1)^2
          data_lm$t3 <- (data_lm$t1)^3
          
          data_lm$t1 <- data_lm$t1/10 # fixed effect estimations reduction in order to better estimate the variances (that are otherwise too small in regards of fixed effects)
          data_lm$t2 <- data_lm$t2/100
          data_lm$t3 <- data_lm$t3/1000
          
        }
        else if(func=="splines"){
          nk = ceiling(length(unique(data_lm$t1))/4)
          noeuds = quantile(data_lm$t1, probs=c(1:(nk))/(nk+1))
          #Bsplines <- as.data.frame(bs(data_lm$t1, knots = noeuds, degree=3, Boundary.knots = range(data_lm$t1), intercept = FALSE))
          NCsplines <- as.data.frame(ns(data_lm$t1, knots = noeuds, Boundary.knots = range(data_lm$t1), intercept = FALSE))
          colnames(NCsplines) <- paste("spline_t",colnames(NCsplines) , sep="")
          NCsplines <- NCsplines*10
          
          data_lm <- cbind.data.frame(data_lm, NCsplines)
          data_lm$t1 <- data_lm$t1/10
          
          splines_DF[gs] <- dim(NCsplines)[2]
          Splines_form <- paste(colnames(NCsplines), collapse=" + ")
        }
        
        if(!separatePatients){
          if(length(levels(data_lm$probe))>1){
          	lmm_H0 <- tryCatch(lmer(expression ~ 1 + (1|probe) + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            if(func=="linear"){
              data_lm$t1 <- data_lm$t1/10
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + (1|probe) + t1 + (t1|probe) + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
          	else if(func=="cubic"){
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + (1|probe) + t1 + t2 + t3 + (t1 + t2 + t3 | probe) + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
                #, control=list(xf.tol=2.2e-14, sing.tol=1e-10))
                #lme(expression ~ 1 + (1|probe) + t1 + t2 + t3, random =  list(as.formula("~ t1 + t2 + t3 | probe"), as.formula("~ 1 | Patient_ID")), method="ML", data=data_lm)
            }
          	else if(func=="splines"){
            	lmm_H1 <- tryCatch(lmer(formula= paste("expression ~ 1 + (1|probe) + ", Splines_form, 
                             " + (", Splines_form, "|probe) + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
          }else{          
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            if(func=="linear"){
              data_lm$t1 <- data_lm$t1/10
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + t1 + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
            else if(func=="cubic"){
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + t1 + t2 + t3 + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
            else if(func=="splines"){
              lmm_H1 <- tryCatch(lmer(formula= paste("expression ~ 1 + ", Splines_form,
                              " + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm, control=list(maxIter=600)),
                       error=function(e){NULL})
            }
          }
        }else{
          if(length(levels(data_lm$probe))>1){
            if(func=="linear"){
              data_lm$t1 <- data_lm$t1/10
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + (1|probe) + t1 + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + (1|probe) + t1 + (t1|Patient_ID) + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
            else if(func=="cubic"){
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + (1|probe) + t1 + t2 + t3 + (1|Patient_ID), REML=FALSE, data=data_lm), #, control=list(xf.tol=2.2e-14, sing.tol=1e-10)),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + (1|probe) + t1 + t2 + t3 + (t1 + t2 + t3 | Patient_ID) + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
            else if(func=="splines"){
              lmm_H0 <- tryCatch(lmer(formula= paste("expression ~ 1 + (1|probe) + ", Splines_form,
                             " + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(formula= paste("expression ~ 1 + (1|probe) + ", Splines_form,
                             " + (", Splines_form, "|Patient_ID) + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
          }else{          
            if(func=="linear"){
              data_lm$t1 <- data_lm$t1/10
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + t1 + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + t1 + (t1|Patient_ID) + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
            else if(func=="cubic"){
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + t1 + t2 + t3 + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + t1 + t2 + t3 + (t1 + t2 + t3 | Patient_ID) + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
            else if(func=="splines"){
              lmm_H0 <- tryCatch(lmer(formula= paste("expression ~ 1 + ", Splines_form,
                             " + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(formula= paste("expression ~ 1 + ", Splines_form,
                             " + (", Splines_form, "|Patient_ID) + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
          }
        }
      }
      else{
        data_temp <- cbind.data.frame(Patient_ID, group.var, TimePoint, expr_temp)
        data_lm <- melt(data_temp, id.vars=c("Patient_ID", "group.var", "TimePoint"), variable.name ="probe", value.name="expression")
        colnames(data_lm)[2] <- "Group"
        colnames(data_lm)[3] <- "t1"
        
        if(func=="cubic"){
          data_lm$t2 <- (data_lm$t1)^2
          data_lm$t3 <- (data_lm$t1)^3
          
          data_lm$t1 <- data_lm$t1/10   # reduction of fixed effect magnitude in order to be able to better estimate variance parameters 
          data_lm$t2 <- data_lm$t2/100  # (otherwise the variances are too small compared to fixed effects)
          data_lm$t3 <- data_lm$t3/1000
          
        }
        else if(func=="splines"){
          nk = ceiling(length(unique(data_lm$t1))/4)
          noeuds = quantile(data_lm$t1, probs=c(1:(nk))/(nk+1))
          #Bsplines <- as.data.frame(bs(data_lm$t1, knots = noeuds, degree=3, Boundary.knots = range(data_lm$t1), intercept = FALSE))
          NCsplines <- as.data.frame(ns(data_lm$t1, knots = noeuds, Boundary.knots = range(data_lm$t1), intercept = FALSE))
          colnames(NCsplines) <- paste("spline_t",colnames(NCsplines) , sep="")
          NCsplines <- NCsplines*10
          
          data_lm <- cbind.data.frame(data_lm, NCsplines)
          data_lm$t1 <- data_lm$t1/10
          
          splines_DF[gs] <- dim(NCsplines)[2]
          Splines_form <- paste(colnames(NCsplines), collapse=" + ")
          SplinesG_form <- paste(paste(colnames(NCsplines), collapse=":Group + "), ":Group", sep="")
        }
        #TODO NA with group.var & splines...
        if(length(levels(data_lm$probe))>1){
          if(func=="linear"){
            data_lm$t1 <- data_lm$t1/10
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + (1|probe) + Group + t1 + (t1|probe) + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + (1|probe) + Group + t1 + t1:Group + (t1:Group|probe) + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
          }
          else if(func=="cubic"){
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + (1|probe) + Group + t1 + t2 + t3 + (t1 + t2 + t3 | probe) + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + (1|probe) + Group + t1 + t1:Group + t2 + t2:Group + t3 + t3:Group + (t1:Group + t2:Group + t3:Group | probe) + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL}) 
          }
          else if(func=="splines"){
          	lmm_H0 <- tryCatch(lmer(formula=paste("expression ~ 1 + (1|probe) + Group + ", Splines_form,
          																				" + (", Splines_form, "|probe) + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm),
          										 error=function(e){NULL})
          	lmm_H1 <- tryCatch(lmer(formula=paste("expression ~ 1 + (1|probe) + Group + ", Splines_form,
          																				" + ", SplinesG_form,
          																				" + (", Splines_form, "|probe)", 
          																				" + (", SplinesG_form, "|probe)",
          																				" + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm),
          										 error=function(e){NULL})
          }
        }else{
          if(func=="linear"){
            data_lm$t1 <- data_lm$t1/10
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + Group + t1 + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + Group + t1 + t1:Group + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
          }
          else if(func=="cubic"){
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + Group + t1 + t2 + t3 + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + Group + t1 + t1:Group + t2 + t2:Group + t3 + t3:Group + (1|Patient_ID), REML=FALSE, data=data_lm),  #, control=list(xf.tol=2.2e-14, sing.tol=1e-10)),
                     error=function(e){NULL})
          }
          else if(func=="splines"){
          	lmm_H0 <- tryCatch(lmer(formula=paste("expression ~ 1 + Group + ", Splines_form,
          																				" + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm, control=list(maxIter=600)),
          										 error=function(e){NULL})
          	lmm_H1 <- tryCatch(lmer(formula=paste("expression ~ 1 + Group + ", Splines_form,
          																				" + ", SplinesG_form, 
          																				" + (1|Patient_ID)", sep=""), REML=FALSE, data=data_lm),
          										 error=function(e){NULL})
          }
        }
      }
		
      if (!is.null(lmm_H0) & !is.null(lmm_H1)) {
        LR[gs] <- lmm_H0@deviance["ML"] - lmm_H1@deviance["ML"]
  	    CVG_H0[gs] <- lmm_H0@dims["cvg"]
  	    CVG_H1[gs] <- lmm_H1@dims["cvg"]
        
        estims <- cbind.data.frame(data_lm, "fitted"=fitted(lmm_H1))
        estims_tab <- acast(data=estims, formula = probe~Patient_ID~t1, value.var="fitted")
        # drop = FALSE by default, which means that missing combination will be kept in the estims_tab and filled with NA
        dimnames(estims_tab)[[3]] <- as.numeric(dimnames(estims_tab)[[3]])*10
        estim_expr[[gs]] <- estims_tab
      } 
      else {
        LR[gs] <- NA
        CVG_H0[gs] <- NA
        CVG_H1[gs] <- NA
        
        estims <- cbind.data.frame(data_lm, "fitted"=NA)
        estims_tab <- acast(data=estims, formula = probe~Patient_ID~t1, value.var="fitted")
        dimnames(estims_tab)[[3]] <- as.numeric(dimnames(estims_tab)[[3]])*10
        estim_expr[[gs]] <- estims_tab
        cat("Unable to fit the mixed models for this gene set\n")
      }
      
      #       "3" = "X-convergence (3)",
      #       "4" = "relative convergence (4)",
      #       "5" = "both X-convergence and relative convergence (5)",
      #       "6" = "absolute function convergence (6)",
      # 
      #       "7" = "singular convergence (7)",
      #       "8" = "false convergence (8)",
      #       "9" = "function evaluation limit reached without convergence (9)",
      #       "10" = "iteration limit reached without convergence (9)",
      #       "14" = "storage has been allocated (?) (14)",
      # 
      #       "15" = "LIV too small (15)",
      #       "16" = "LV too small (16)",
      #       "63" = "fn cannot be computed at initial par (63)",
      #       "65" = "gr cannot be computed at initial par (65)")
    
    }
    else{
	    LR[gs] <- NA
	    CVG_H0[gs] <- NA
	    CVG_H1[gs] <- NA
	    
	    estims <- cbind.data.frame(data_lm, "fitted"=NA)
	    estims_tab <- acast(data=estims, formula = probe~Patient_ID~t1, value.var="fitted")
	    dimnames(estims_tab)[[3]] <- as.numeric(dimnames(estims_tab)[[3]])*10
	    estim_expr[[gs]] <- estims_tab
	    cat("The size of the gene set is problematic (too many or too few genes)\n")
	}
    
    cat(paste(gs,"/", length(gmt$genesets)," gene sets analyzed\n", sep=""))
  }
  tcgsa <- list("fit"=as.data.frame(cbind(LR, CVG_H0, CVG_H1)), "func_form"=func, "GeneSets_gmt"=gmt, 
  							"group.var"=group.var, "separatePatients"=separatePatients, "Estimations"=estim_expr, 
  							"splines_DF"=ifelse(length(which(!is.na(unique(splines_DF))))>0, max(unique(splines_DF), na.rm=T), NA)
  							)
  class(tcgsa) <- "TcGSA"
  return(tcgsa)
}
