plot.ClusteredTrends <- function(x, ...){
  maxK <- x$MaxNbClust
  f <- factor(x$NbClust, levels=c(1:maxK))
  barplot(height=summary(f),
          xlab="Number of distinct trends", ylab= "Number of gene sets",
          col=rainbow(maxK),
  				main=paste(formatC(mean(x$NbClust), digits=3), "trends by significant gene set (on average)"),
  				ylim=c(0, sum(summary(f)))  				
  )      
}