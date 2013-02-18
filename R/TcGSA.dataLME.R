#'@keywords internal

TcGSA.dataLME<-
	function(expr, design, subject_name="Patient_ID", time_name="TimePoint", 
					 covariates_fixed="", time_covariates="",
					 group_name=""){

	data_temp <- cbind.data.frame(design, expr)
	
	ID_vars <- c(subject_name, time_name, covariates_fixed, time_covariates, group_name)
	if(length(which(ID_vars==""))>0){ID_vars <- ID_vars[-which(ID_vars=="")]}
	
	data_lm <- melt(data_temp, id.vars=ID_vars, variable.name ="probe", value.name="expression", measure.vars=colnames(expr))
	colnames(data_lm)[which(colnames(data_lm)==time_name)] <- "t1"

	data_lm$t1 <- data_lm$t1/10 # fixed effect estimations reduction in order to better estimate the variances (that are otherwise too small in regards of fixed effects)
	data_lm$t2 <- (data_lm$t1)^2
	data_lm$t3 <- (data_lm$t1)^3

	nk = ceiling(length(unique(design[,time_name]))/4)
	noeuds = quantile(design[,time_name], probs=c(1:(nk))/(nk+1))
	NCsplines <- as.data.frame(ns(design[,time_name], knots = noeuds, Boundary.knots = range(design[,time_name]), intercept = FALSE))
	colnames(NCsplines) <- paste("spline_t",colnames(NCsplines) , sep="")
	NCsplines <- NCsplines*10
	data_lm <- cbind.data.frame(data_lm, NCsplines)
	
	return(data_lm)
}
	