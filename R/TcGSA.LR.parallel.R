TcGSA.LR.parallel <-
function(Nproc, type_connec, expr, gmt, Patient_ID, TimePoint, func = "linear", maxGSsize=500, group.var=NULL, separatePatients=FALSE, monitorfile=""){
	
#   library(GSA)
#   library(lme4)
#   library(reshape2)
#   library(splines) 
   library(doSNOW)
	

  if(!is.null(group.var) & separatePatients){
    stop("'separatePatients' is TRUE while 'group.var' is not NULL.\n This is an attempt to separate patients in a multiple group setting.\n This is not handled by the TcGSA.LR function.\n\n")
  }
  
   
  cl <- makeCluster(Nproc, type = type_connec)
  registerDoSNOW(cl)
  
  res_par <- foreach(gs=1:length(gmt$genesets), .packages=c("lme4", "reshape2", "splines")) %dopar% {
    probes <- intersect(gmt$genesets[[gs]], rownames(expr))

    
    if(length(probes)>0 & length(probes)<maxGSsize){                                                       
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
          
        }else if(func=="splines"){
          nk = ceiling(length(unique(data_lm$t1))/4)
          noeuds = quantile(data_lm$t1, probs=(c(0:(nk+1))/(nk+1))[-c(1,(nk+1+1))])
          Bsplines <- as.data.frame(bs(data_lm$t1, knots = noeuds, Boundary.knots = range(data_lm$t1), intercept = FALSE))
          colnames(Bsplines) <- paste("spline_t",colnames(Bsplines) , sep="")
          Bsplines <- Bsplines*10
          data_lm <- cbind.data.frame(data_lm, Bsplines)
        }
        
        
        
        if(!separatePatients){
          if(length(levels(data_lm$probe))>1){
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + probe + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            if(func=="linear"){
              data_lm$t1 <- data_lm$t1/10
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + probe + t1 + (t1|probe) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }else if(func=="cubic"){
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + probe + t1 + t2 + t3 + (t1 + t2 + t3 | probe) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }else if(func=="splines"){
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + probe + spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5
                             + (spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5|probe) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
          }else{          
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            if(func=="linear"){
              data_lm$t1 <- data_lm$t1/10
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + t1 + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }else if(func=="cubic"){
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + t1 + t2 + t3 + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }else if(func=="splines"){
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5
                              + (1|Patient_ID), REML=FALSE, data=data_lm, control=list(maxIter=600)),
                       error=function(e){NULL})
            }
          }
        }else{
          if(length(levels(data_lm$probe))>1){
            if(func=="linear"){
              data_lm$t1 <- data_lm$t1/10
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + probe + t1 + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + probe + t1 + (t1|Patient_ID) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }else if(func=="cubic"){
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + probe + t1 + t2 + t3 + (1|Patient_ID:probe), REML=FALSE, data=data_lm), #, control=list(xf.tol=2.2e-14, sing.tol=1e-10)),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + probe + t1 + t2 + t3 + (t1 + t2 + t3 | Patient_ID) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),  #, control=list(xf.tol=2.2e-14, sing.tol=1e-10)),
                       error=function(e){NULL})
            }else if(func=="splines"){
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + probe + spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5
                             + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + probe + spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5
                             + (spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5|Patient_ID) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
          }else{          
            if(func=="linear"){
              data_lm$t1 <- data_lm$t1/10
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + t1 + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + t1 + (t1|Patient_ID) + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }else if(func=="cubic"){
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + t1 + t2 + t3 + (1|Patient_ID), REML=FALSE, data=data_lm), #, control=list(xf.tol=2.2e-14, sing.tol=1e-10)),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + t1 + t2 + t3 + (t1 + t2 + t3 | Patient_ID) + (1|Patient_ID), REML=FALSE, data=data_lm), #, control=list(xf.tol=2.2e-14, sing.tol=1e-10)),
                       error=function(e){NULL})
            }else if(func=="splines"){
              lmm_H0 <- tryCatch(lmer(expression ~ 1 + spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5
                             + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
              lmm_H1 <- tryCatch(lmer(expression ~ 1 + spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5
                             + (spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5|Patient_ID) + (1|Patient_ID), REML=FALSE, data=data_lm),
                       error=function(e){NULL})
            }
          }
        }
      }else{
        data_temp <- cbind.data.frame(Patient_ID, group.var, TimePoint, expr_temp)
        data_lm <- melt(data_temp, id.vars=c("Patient_ID", "group.var", "TimePoint"), variable.name ="probe", value.name="expression")
        colnames(data_lm)[2] <- "Group"
        colnames(data_lm)[3] <- "t1"
        
        if(func=="cubic"){
          data_lm$t2 <- (data_lm$t1)^2
          data_lm$t3 <- (data_lm$t1)^3
          
          data_lm$t1 <- data_lm$t1/10 #reduction des effets fixes afin de pouvoir mieux estimer les variances (sinon variances trop faibles aux regards des effets fixes)
          data_lm$t2 <- data_lm$t2/100
          data_lm$t3 <- data_lm$t3/1000
          
        }else if(func=="splines"){
          nk = ceiling(length(unique(data_lm$t1))/4)
          noeuds = quantile(data_lm$t1, probs=(c(0:(nk+1))/(nk+1))[-c(1,(nk+1+1))])
          Bsplines <- as.data.frame(bs(data_lm$t1, knots = noeuds, Boundary.knots = range(data_lm$t1), intercept = FALSE))
          colnames(Bsplines) <- paste("spline_t",colnames(Bsplines) , sep="")
          Bsplines <- Bsplines*10
          data_lm <- cbind.data.frame(data_lm, Bsplines)
        }
        
        if(length(levels(data_lm$probe))>1){
          if(func=="linear"){
            data_lm$t1 <- data_lm$t1/10
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + probe + Group + t1 + (t1|probe) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + probe + Group + t1 + t1:Group + (t1:Group|probe) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
          }else if(func=="cubic"){
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + probe + Group + t1 + t2 + t3 + (t1 + t2 + t3 | probe) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                     error=function(e){NULL}) 
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + probe + Group + t1 + t1:Group + t2 + t2:Group + t3 + t3:Group + (t1:Group + t2:Group + t3:Group | probe) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                     error=function(e){NULL}) 
          }else if(func=="splines"){
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + probe + Group + spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5
                           + (spline_t1 + spline_t2 + spline_t3 + spline_t4 +spline_t5|probe) + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + probe + Group + spline_t1 + spline_t2 + spline_t3 + spline_t4 + spline_t5
                           + spline_t1:Group + spline_t2:Group + spline_t3:Group + spline_t4:Group +spline_t5:Group
                           + (spline_t1 + spline_t2 + spline_t3 + spline_t4 + spline_t5|probe) 
                           + (spline_t1:Group + spline_t2:Group + spline_t3:Group + spline_t4:Group +spline_t5:Group|probe)
                           + (1|Patient_ID:probe), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
          }
        }else{
          if(func=="linear"){
            data_lm$t1 <- data_lm$t1/10
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + Group + t1 + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + Group + t1 + t1:Group + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
          }else if(func=="cubic"){
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + Group + t1 + t2 + t3 + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + Group + t1 + t1:Group + t2 + t2:Group + t3 + t3:Group + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL}) 
          }else if(func=="splines"){
            lmm_H0 <- tryCatch(lmer(expression ~ 1 + Group + spline_t1 + spline_t2 + spline_t3 + spline_t4 + spline_t5
                           + (1|Patient_ID), REML=FALSE, data=data_lm, control=list(maxIter=600)),
                     error=function(e){NULL})
            lmm_H1 <- tryCatch(lmer(expression ~ 1 + Group + spline_t1 + spline_t2 + spline_t3 + spline_t4 + spline_t5
                           + spline_t1:Group + spline_t2:Group + spline_t3:Group + spline_t4:Group +spline_t5:Group
                           + (1|Patient_ID), REML=FALSE, data=data_lm),
                     error=function(e){NULL})
          }
        }
      }
	  #lmm_fits_H0[[gs]] <- lmm_H0
	  #lmm_fits_H1[[gs]] <- lmm_H1
		
      if (!is.null(lmm_H0) & !is.null(lmm_H1)) {
      	LR <- lmm_H0@deviance["ML"] - lmm_H1@deviance["ML"]
      	CVG_H0 <- lmm_H0@dims["cvg"]
      	CVG_H1 <- lmm_H1@dims["cvg"]
    
      	FixEf <- fixef(lmm_H1)
      	RanEf <- ranef(lmm_H1)
      } else {
      	LR <- NA
      	CVG_H0 <- NA
      	CVG_H1 <- NA
      	
      	FixEf <- NA
      	RanEf <- NA
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
      #       LR[gs] <- lmm_H0@deviance["ML"] - lmm_H1@deviance["ML"]
    
    }else{
	    LR <- NA
	    CVG_H0 <- NA
	    CVG_H1 <- NA
	    FixEf <- NA
	    RanEf <- NA
	}
    line_number <- 0
    try(line_number <- length(readLines(monitorfile)), silent=TRUE)
    cat(paste(line_number+1,"/", length(gmt$genesets)," gene sets analyzed (geneset ", gs, ")\n", sep=""), file=monitorfile, append = TRUE)

    res <- list("LR"=LR, "CVG_H0"=CVG_H0, "CVG_H1"=CVG_H1, "FixEf"=FixEf, "RanEf"=RanEf)
  }
  
  LR <- numeric(length(gmt$genesets))
  CVG_H0 <- numeric(length(gmt$genesets))
  CVG_H1 <- numeric(length(gmt$genesets))
  FixEf <- list()
  RanEf <- list()
  
  cat("Combining the results...")
  for (gs in 1:length(gmt$genesets)){
    LR[gs] <- res_par[[gs]][["LR"]]
    CVG_H0[gs] <- res_par[[gs]][["CVG_H0"]]
    CVG_H1[gs] <- res_par[[gs]][["CVG_H1"]]
    FixEf[[gs]] <- res_par[[gs]][["FixEf"]]
    RanEf[[gs]] <- res_par[[gs]][["RanEf"]]
  }
  tcgsa <- list("fit"=as.data.frame(cbind(LR, CVG_H0, CVG_H1)), "func_form"=func, "GeneSets_gmt"=gmt, "group.var"=group.var, "separatePatients"=separatePatients, "Estimations"=list("FixEf"=FixEf, "RanEf"=RanEf))
  class(tcgsa) <- "TcGSA"
  stopCluster(cl)
  return(tcgsa)
}
