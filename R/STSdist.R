#'@keywords internalSTSdist <- function(m, time){

STSdist <- function(m, time){
	nr <- dim(m)[1]
	nc <- dim(m)[2]
	mm <- matrix(NA, nrow=nr, ncol=nr)
	mdiff <- (m[, -1] - m[,-nc])/(time[-1] - time[-nc])
	for(i in 1:(nr-1)){
		for(j in (i+1):nr){
			mm[j,i] <- sqrt(sum((mdiff[i,] -mdiff[j,])^2))
		}
	}
	return(as.dist(mm))
}

