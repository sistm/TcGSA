print.ClusteredTrends <- function(x, ...){
  maxK <- x$MaxNbClust
  f <- factor(x$NbClust, levels=c(1:maxK))
  levels(f)[-1] <- paste(levels(f)[-1], "trends")
  levels(f)[1] <- paste(levels(f)[1], "trend")
  cat("\t\t\tA ClusteredTrends object")
  cat("\n")
  cat("\n")
  cat("Distribution of the number of trends per gene sets:")
  cat("\n\t")
  for(i in 1:maxK){
    cat(levels(f)[i])
    cat(": ")
    if(i==1){cat(" ")}
    cat(summary(f)[i])
    cat("\n\t")
  }
  cat("Total number of trends:", sum(x$NbClust), "(out of", length(x$NbClust), "significant gene sets)\n") 
  cat("\n")
  
  cat("Maximal number of clusters tested:", maxK, "\n")

  cat("\n")
  cat("Mean number of trends by significant gene set:", formatC(mean(x$NbClust), digits=3), "\n") 
}