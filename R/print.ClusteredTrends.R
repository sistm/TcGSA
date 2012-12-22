print.ClusteredTrends <- function(x, ...){
  maxK <- x$MaxNbClust
  f <- factor(x$NbClust, levels=c(1:maxK))
  levels(f)[-1] <- paste(levels(f)[-1], "clusters")
  levels(f)[1] <- paste(levels(f)[1], "cluster")
  cat("\t\t\tA ClusteredTrends object")
  cat("\n")
  cat("\n")
  cat("Distribution of the number of clusters per gene sets:")
  cat("\n\t")
  for(i in 1:maxK){
    cat(levels(f)[i])
    cat(": ")
    if(i==1){cat(" ")}
    cat(summary(f)[i])
    cat("\n\t")
  }
  cat("\n")
  cat("Maximal number of clusters tested:")
  cat("\n\t")
  cat(maxK)
}