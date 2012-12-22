plot.ClusteredTrends <- function(x, ...){
  maxK <- x$MaxNbClust
  f <- factor(x$NbClust, levels=c(1:maxK))
  barplot(height=summary(f),
          xlab="Number of clusters", ylab= "Number of gene sets",
          col=rainbow(maxK)
  )      
}