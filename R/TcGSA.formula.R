#'@keywords internal
TcGSA.formula <-
	function(design, subject_name="Patient_ID", time_name="TimePoint",
			 covariates_fixed="", time_covariates="", group_name="",
			 separateSubjects=FALSE, crossedRandom_fixed=FALSE, time_crossedRandom=FALSE,
			 time_func = "linear"){
		# TIME function
		if(time_func=="linear"){
			time <- "t1"
		}else if(time_func=="cubic"){
			time <- "t1 + t2 + t3"
		}else if(time_func=="splines"){
			nk = ceiling(length(unique(design[,time_name]))/4)
			noeuds = quantile(design[,time_name], probs=c(1:(nk))/(nk+1))
			NCsplines <- as.data.frame(ns(design[,time_name], knots = noeuds, Boundary.knots = range(design[,time_name]), intercept = FALSE))
			time <- paste(" + spline_t", colnames(NCsplines), collapse="", sep="")
			#time <- paste(" + spline_t",1:(nk+1) , sep="", collapse="")
		}else{
			time <- time_func
			#a user specified function of time.
			#This must be an expression involving only columnnames of the design matrix
		}
		if(time_func %in% c("linear", "cubic", "splines")){
			time_split <- str_trim(str_split(paste("+", time, collapse=" "), "\\+")[[1]])
			if(length(which(time_split==""))>0){time_split <- time_split[-which(time_split=="")]}
			time_DF <- length(time_split)
		}else if((time_func %in% colnames(design)) && is.factor(design[, time_func])){
			time_DF <- length(levels(design[, time_func]))-1
		}else {
			time_split <- gsub(" ", "", unlist(lapply(unlist(lapply(unlist(str_split(time_func, "\\+")), FUN=str_split, pattern="\\/")), FUN=str_split, pattern="\\*")))
			time_DF <- 0
			for (v in time_split){
				time_DF <- time_DF + ifelse(is.numeric(design[, v]), 1, length(levels(as.factor(design[, v])))-1)
			}
		}
		# Covariates
		if(covariates_fixed[1]!=""){
			covariates_fixed <- paste(" + ", covariates_fixed, collapse="", sep="")
		}
		if(time_covariates[1]!=""){
			tc <- NULL
			for (c in time_covariates){
				tc <- paste(tc, paste(" +", paste(time_split, c, sep=":", collapse=" + "), collapse=" "), sep="")
			}
			time_covariates <- tc
		}
		if(time_func!="cubic"){
			if (!crossedRandom_fixed){
				if(!time_crossedRandom){
					if(group_name=="" & !separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", 
										   covariates_fixed, time_covariates, sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ",time, time_covariates,
										   " + (0 + ", time, "|probe) ", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ")", sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ",time, time_covariates,
												  " + (1|", subject_name, ")", sep="")
					}else if(group_name=="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe)", covariates_fixed,
										   " + (1|", subject_name, ") + ", time, time_covariates, " + (0 + ", time, "|probe)",sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe)+ (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", time, "|probe) + (0 + ", time, "|", subject_name, ")", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", time, "|", subject_name, ") + (1|", subject_name, ")", sep="")
					}else if(group_name!="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", time, "|probe) + ", time, ":", group_name, sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", time, "|probe) + (0 + ", time, "|", subject_name, ") + ", time, ":", group_name, sep="")
						formula_H0_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + ", time, ":", group_name, sep="")
						formula_H1_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", time, "|", subject_name, ") + ", time, ":", group_name, sep="")
					}
				}else{
					if(group_name=="" & !separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe) ",covariates_fixed, 
										   " + (1|", subject_name, ")", sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe)", covariates_fixed, " + ",time,time_covariates,
										   " + (0 + ", time, "|probe) + (1|", subject_name, ")", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ")", sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ",time, time_covariates,
												  " + (1|", subject_name, ")", sep="")
					}else if(group_name=="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe)", covariates_fixed,
										   " + (1|", subject_name, ") + ", time, time_covariates, " + ",time,":probe",sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe)+ (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + ",time,":probe + (0 + ", time, "|",subject_name,":probe)", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", time, "|", subject_name, ") + (1|", subject_name, ")", sep="")
					}else if(group_name!="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + ", time, ":probe + ", time, ":", group_name, sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + ", time, ":probe + (0 + ", time, "|probe:", subject_name, ") + ", time, ":", group_name, sep="")
						formula_H0_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + ", time, ":", group_name, sep="")
						formula_H1_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", time, "|", subject_name, ") + ", time, ":", group_name, sep="")
					}
				}
			}else{
				if(!time_crossedRandom){
					if(group_name=="" & !separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe", covariates_fixed, time_covariates,
										   " + (1|", subject_name, ":probe)", sep="")
						formula_H1 = paste("expression ~ 1 + probe", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", time, "|probe) + (1|", subject_name, ":probe)", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ")", sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time,time_covariates,
												  " + (1|", subject_name, ")", sep="")
					}else if(group_name=="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe", covariates_fixed,
										   " + (1|", subject_name, ":probe) + ", time, time_covariates, " + (0 + ", time, "|probe)", sep="")
						formula_H1 = paste("expression ~ 1 + probe", covariates_fixed,
										   " + (1|", subject_name, ":probe) + ", time, time_covariates, " + (0 + ", time, "|probe) + (0 + ", time, "|", subject_name, ")", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", time, "|", subject_name, ") + (1|", subject_name, ")", sep="")
					}else if(group_name!="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe + (1|", subject_name, ":probe)", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", time, "|probe) + ", time, ":", group_name, sep="")
						formula_H1 = paste("expression ~ 1 + probe + (1|", subject_name, ":probe)", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", time, "|probe) + (0 + ", substr(time,1,2), "|", subject_name,") + (0 + ", substr(time,6,7), "|", subject_name,") + (0 + ", substr(time,11,12), "|", subject_name,") + ", time, ":", group_name, sep="")
						formula_H0_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + ", time, ":", group_name, sep="")
						formula_H1_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", substr(time,1,2), "|", subject_name,") + (0 + ", substr(time,6,7), "|", subject_name,") + (0 + ", substr(time,11,12), "|", subject_name,") + ", time, ":", group_name, sep="")
					}
				}else{
					if(group_name=="" & !separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe", covariates_fixed,time_covariates,
										   " + (1|", subject_name, ":probe)", sep="")
						formula_H1 = paste("expression ~ 1 + probe", covariates_fixed," + ", time,time_covariates,
										   " + (0 + ", time, "|probe) + (1|", subject_name, ":probe)", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ")", sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time,time_covariates,
												  " + (1|", subject_name, ")", sep="")
					}else if(group_name=="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe", covariates_fixed,
										   " + (1|", subject_name, ":probe) + ", time, time_covariates, " + " ,time, ":probe", sep="")
						formula_H1 = paste("expression ~ 1 + probe", covariates_fixed,
										   " + (1|", subject_name, ":probe) + ", time, time_covariates, " + " ,time, ":probe + (0 + ", time, "|", subject_name,":probe)", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, " + (0 + ", substr(time,1,2), "|", subject_name,") + (0 + ", substr(time,6,7), "|", subject_name,") + (0 + ", substr(time,11,12), "|", subject_name,")", sep="")
					}else if(group_name!="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe + (1|", subject_name, ":probe)", covariates_fixed, " + ", time, time_covariates,
										   " + ", time, ":probe + ", time, ":", group_name, sep="")
						formula_H1 = paste("expression ~ 1 + probe + (1|", subject_name, ":probe)", covariates_fixed, " + ", time, time_covariates,
										   " + ", time, ":probe + (0 + ", time, "|probe:", subject_name, ") + ", time, ":", group_name, sep="")
						formula_H0_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + ", time, ":", group_name, sep="")
						formula_H1_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", time, "|", subject_name, ") + ", time, ":", group_name, sep="")
					}
				}
			}
		}else{
			if (!crossedRandom_fixed){
				if(!time_crossedRandom){
					if(group_name=="" & !separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", 
										   covariates_fixed, time_covariates, sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ",time, time_covariates,
										   " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  ", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ")", sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ",time, time_covariates,
												  " + (1|", subject_name, ")", sep="")
					}else if(group_name=="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe)", covariates_fixed,
										   " + (1|", subject_name, ") + ", time, time_covariates, " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe) ",sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe)+ (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  + (0 + ", substr(time,1,2), "|", subject_name, ") + (0 + ", substr(time,6,7), "|", subject_name, ") + (0 + ", substr(time,11,12), "|", subject_name, ")", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", substr(time,1,2), "|",subject_name,") + (0 + ", substr(time,6,7), "|",subject_name,") + (0 + ", substr(time,11,12), "|",subject_name,") + (1|", subject_name, ")", sep="")
					}else if(group_name!="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  + (0 + ", substr(time,1,2), "|",subject_name,") + (0 + ", substr(time,6,7), "|",subject_name,") + (0 + ", substr(time,11,12), "|",subject_name,") + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H0_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H1_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", substr(time,1,2), "|",subject_name,") + (0 + ", substr(time,6,7), "|",subject_name,") + (0 + ", substr(time,11,12), "|",subject_name,") + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
					}
				}else{
					if(group_name=="" & !separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe) ",covariates_fixed, 
										   " + (1|", subject_name, ")", sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe)", covariates_fixed, " + ",time,time_covariates,
										   " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  + (1|", subject_name, ")", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ")", sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ",time, time_covariates,
												  " + (1|", subject_name, ")", sep="")
					}else if(group_name=="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe)", covariates_fixed,
										   " + (1|", subject_name, ") + ", time, time_covariates, " + ",substr(time,1,2),":probe + ", substr(time,6,7),":probe + ", substr(time,11,12),":probe",sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe)+ (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + ",substr(time,1,2),":probe + ", substr(time,6,7),":probe + ", substr(time,11,12),":probe + (0 + ", substr(time,1,2), "|",subject_name,":probe) + (0 + ", substr(time,6,7), "|",subject_name,":probe) + (0 + ", substr(time,11,12), "|",subject_name,":probe)", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", substr(time,1,2), "|",subject_name,") + (0 + ", substr(time,6,7), "|",subject_name,") + (0 + ", substr(time,11,12), "|",subject_name,") + (1|", subject_name, ")", sep="")
					}else if(group_name!="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + ", substr(time,1,2), ":probe + ", substr(time,6,7), ":probe + ", substr(time,11,12), ":probe + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H1 = paste("expression ~ 1 + (1|probe) + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
										   " + ", substr(time,1,2), ":probe + ", substr(time,6,7), ":probe + ", substr(time,11,12), ":probe + (0 + ", substr(time,1,2), "|probe:", subject_name, ") + (0 + ", substr(time,6,7), "|probe:", subject_name, ") + (0 + ", substr(time,11,12), "|probe:", subject_name, ") + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H0_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H1_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", substr(time,1,2), "|",subject_name,") + (0 + ", substr(time,6,7), "|",subject_name,") + (0 + ", substr(time,11,12), "|",subject_name,") + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
					}
				}
			}else{
				if(!time_crossedRandom){
					if(group_name=="" & !separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe", covariates_fixed, time_covariates,
										   " + (1|", subject_name, ":probe)", sep="")
						formula_H1 = paste("expression ~ 1 + probe", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  + (1|", subject_name, ":probe)", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ")", sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time,time_covariates,
												  " + (1|", subject_name, ")", sep="")
					}else if(group_name=="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe", covariates_fixed,
										   " + (1|", subject_name, ":probe) + ", time, time_covariates, " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe) ", sep="")
						formula_H1 = paste("expression ~ 1 + probe", covariates_fixed,
										   " + (1|", subject_name, ":probe) + ", time, time_covariates, " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  + (0 + ", substr(time,1,2), "|", subject_name, ") + (0 + ", substr(time,6,7), "|", subject_name, ") + (0 + ", substr(time,11,12), "|", subject_name, ")", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", substr(time,1,2), "|",subject_name,") + (0 + ", substr(time,6,7), "|",subject_name,") + (0 + ", substr(time,11,12), "|",subject_name,") + (1|", subject_name, ")", sep="")
					}else if(group_name!="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe + (1|", subject_name, ":probe)", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H1 = paste("expression ~ 1 + probe + (1|", subject_name, ":probe)", covariates_fixed, " + ", time, time_covariates,
										   " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  + (0 + ", substr(time,1,2), "|", subject_name,") + (0 + ", substr(time,6,7), "|", subject_name,") + (0 + ", substr(time,11,12), "|", subject_name,") + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H0_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H1_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", substr(time,1,2), "|", subject_name,") + (0 + ", substr(time,6,7), "|", subject_name,") + (0 + ", substr(time,11,12), "|", subject_name,") + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
					}
				}else{
					if(group_name=="" & !separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe", covariates_fixed,time_covariates,
										   " + (1|", subject_name, ":probe)", sep="")
						formula_H1 = paste("expression ~ 1 + probe", covariates_fixed," + ", time,time_covariates,
										   " + (0 + ", substr(time,1,2), "|probe) + (0 + ", substr(time,6,7), "|probe) + (0 + ", substr(time,11,12), "|probe)  + (1|", subject_name, ":probe)", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ")", sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed, " + ", time,time_covariates,
												  " + (1|", subject_name, ")", sep="")
					}else if(group_name=="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe", covariates_fixed,
										   " + (1|", subject_name, ":probe) + ", time, time_covariates, " + " ,substr(time,1,2), ":probe + ", substr(time,6,7), ":probe + ", substr(time,11,12), ":probe", sep="")
						formula_H1 = paste("expression ~ 1 + probe", covariates_fixed,
										   " + (1|", subject_name, ":probe) + ", time, time_covariates, " + " ,substr(time,1,2), ":probe + ", substr(time,6,7), ":probe + ", substr(time,11,12), ":probe + (0 + ", substr(time,1,2), "|", subject_name,":probe) + (0 + ", substr(time,6,7), "|", subject_name,":probe) + (0 + ", substr(time,11,12), "|", subject_name,":probe)", sep="")
						formula_H0_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, sep="")
						formula_H1_1probe = paste("expression ~ 1", covariates_fixed,
												  " + (1|", subject_name, ") + ", time, time_covariates, " + (0 + ", substr(time,1,2), "|", subject_name,") + (0 + ", substr(time,6,7), "|", subject_name,") + (0 + ", substr(time,11,12), "|", subject_name,")", sep="")
					}else if(group_name!="" & separateSubjects){
						formula_H0 = paste("expression ~ 1 + probe + (1|", subject_name, ":probe)", covariates_fixed, " + ", time, time_covariates,
										   " + ", substr(time,1,2), ":probe + ", substr(time,6,7), ":probe + ", substr(time,11,12), ":probe + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H1 = paste("expression ~ 1 + probe + (1|", subject_name, ":probe)", covariates_fixed, " + ", time, time_covariates,
										   " + ", substr(time,1,2), ":probe + ", substr(time,6,7), ":probe + ", substr(time,11,12), ":probe + (0 + ", substr(time,1,2), "|probe:", subject_name, ") + (0 + ", substr(time,6,7), "|probe:", subject_name, ") + (0 + ", substr(time,11,12), "|probe:", subject_name, ") + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H0_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
						formula_H1_1probe = paste("expression ~ 1 + (1|", subject_name, ")", covariates_fixed, " + ", time, time_covariates,
												  " + (0 + ", substr(time,1,2), "|",subject_name,") + (0 + ", substr(time,6,7), "|",subject_name,") + (0 + ", substr(time,11,12), "|",subject_name,") + ", substr(time,1,2), ":", group_name, " + ", substr(time,6,7), ":", group_name, " + ", substr(time,11,12), ":", group_name, sep="")
					}
				}
			}
		}
		return(list("H0"=c("reg"=formula_H0, "1probe"=formula_H0_1probe), "H1"=c("reg"=formula_H1, "1probe"=formula_H1_1probe), "time_DF"=time_DF))	
	}
