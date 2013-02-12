TcGSA.LR.new <-
function(expr, gmt, design, subject_name="Patient_ID", time_name="TimePoint", crossedRandom=FALSE,
				 covariates_fixed="", time_covariates="",
				 time_func = "linear", group.var=NULL, separateSubjects=FALSE,
				 minGSsize=10, maxGSsize=500){
  
  #DEBUG: data(data_simu_TcGSA);
  #expr=expr_1grp; Patient_ID=design$Patient_ID; TimePoint=design$TimePoint; gmt=gmt_sim; time_func="linear"; gs=2; maxGSsize=500; group.var=NULL;
  #for (i in which(tcGSA_cubic_modules$CVG_H1>6)){
  #  print(paste(i, ":", length(intersect(gmt$genesets[[i]], rownames(expr)))))
  #}
  
  if(!is.null(group.var) && separateSubjects){
    stop("'separateSubjects' is TRUE while 'group.var' is not NULL.\n This is an attempt to separate subjects in a multiple group setting.\n This is not handled by the TcGSA.LR function.\n\n")
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
  
  my_formul <- TcGSA.formula(design=design, subject_name=subject_name, time_name=time_name, crossedRandom=crossedRandom, 
  													 covariates_fixed=covariates_fixed, time_covariates=time_covariates,
  													 time_func=time_func, group.var=group.var, separateSubjects=separateSubjects)
  
  for (gs in 1:length(gmt$genesets)){
    probes <- intersect(gmt$genesets[[gs]], rownames(expr))
    
    if(length(probes)>0 && length(probes)<=maxGSsize && length(probes)>=minGSsize){                                                       
    	expr_temp <- t(expr[probes, ])
    	rownames(expr_temp) <- NULL
    	data_lme  <- TcGSA.dataLME(expr=expr_temp, design=design, subject_name=subject_name, time_name=time_name, 
    														 covariates_fixed=covariates_fixed, time_covariates=time_covariates,
    														 group.var=group.var)
    	
      if(length(levels(data_lme$probe))>1){
          lmm_H0 <- tryCatch(lmer(formula =my_formul[["H0"]]["reg"], REML=FALSE, data=data_lme),
                   error=function(e){NULL})
          lmm_H1 <- tryCatch(lmer(formula =my_formul[["H1"]]["reg"], REML=FALSE, data=data_lme),
                   error=function(e){NULL})
      }
    	else{
    		lmm_H0 <- tryCatch(lmer(formula =my_formul[["H0"]]["1probe"], REML=FALSE, data=data_lme),
    											 error=function(e){NULL})
    		lmm_H1 <- tryCatch(lmer(formula =my_formul[["H1"]]["1probe"], REML=FALSE, data=data_lme),
    											 error=function(e){NULL})
      }
		
      if (!is.null(lmm_H0) & !is.null(lmm_H1)) {
        LR[gs] <- lmm_H0@deviance["ML"] - lmm_H1@deviance["ML"]
  	    CVG_H0[gs] <- lmm_H0@dims["cvg"]
  	    CVG_H1[gs] <- lmm_H1@dims["cvg"]
        
        estims <- cbind.data.frame(data_lme, "fitted"=fitted(lmm_H1))
        estims_tab <- acast(data=estims, formula = probe~Patient_ID~t1, value.var="fitted")
        # drop = FALSE by default, which means that missing combination will be kept in the estims_tab and filled with NA
        dimnames(estims_tab)[[3]] <- as.numeric(dimnames(estims_tab)[[3]])*10
        estim_expr[[gs]] <- estims_tab
      } 
      else {
        LR[gs] <- NA
        CVG_H0[gs] <- NA
        CVG_H1[gs] <- NA
        
        estims <- cbind.data.frame(data_lme, "fitted"=NA)
        estims_tab <- acast(data=estims, formula = probe~Patient_ID~t1, value.var="fitted")
        dimnames(estims_tab)[[3]] <- as.numeric(dimnames(estims_tab)[[3]])*10
        estim_expr[[gs]] <- estims_tab
        cat("Unable to fit the mixed models for this gene set\n")
      }
      
    	#		CONVERGENCE DIAGNOSTICS IN LME4
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
    	#
    	
    }
    else{
	    LR[gs] <- NA
	    CVG_H0[gs] <- NA
	    CVG_H1[gs] <- NA
	    
	    estims <- cbind.data.frame(data_lme, "fitted"=NA)
	    estims_tab <- acast(data=estims, formula = probe~Patient_ID~t1, value.var="fitted")
	    dimnames(estims_tab)[[3]] <- as.numeric(dimnames(estims_tab)[[3]])*10
	    estim_expr[[gs]] <- estims_tab
	    cat("The size of the gene set is problematic (too many or too few genes)\n")
	}
    
    cat(paste(gs,"/", length(gmt$genesets)," gene sets analyzed\n", sep=""))
  }
  tcgsa <- list("fit"=as.data.frame(cbind(LR, CVG_H0, CVG_H1)), "func_form"=time_func, "GeneSets_gmt"=gmt, 
  							"group.var"=group.var, "separatePatients"=separateSubjects, "Estimations"=estim_expr, 
  							"splines_DF"=ifelse(length(which(!is.na(unique(splines_DF))))>0, max(unique(splines_DF), na.rm=T), NA)
  							)
  class(tcgsa) <- "TcGSA"
  return(tcgsa)
}
