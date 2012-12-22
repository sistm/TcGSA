plot1GS <- 
  function(expr, gmt, Patient_ID,TimePoint, geneset.name, 
            baseline=NULL,
            group.var=NULL, Group_ID_paired=NULL, ref=NULL, group_of_interest=NULL,
            FUNcluster=NULL, clustering_metric="euclidian", clustering_method="ward", B=500,
            max_trends=4, aggreg.fun="median", trend.fun="median",
            methodOptiClust = "firstSEmax",
            indiv="genes",
            verbose=TRUE,
            clustering=TRUE, showTrend=TRUE, smooth=TRUE,
            time_unit="", title=NULL, desc=TRUE,
            lab.cex=1, axis.cex=1, main.cex=1, y.lab.angle=90, x.axis.angle=45,
            y.lim=NULL, x.lim=NULL, 
            gg.add=theme()
           ){
  
  library(ggplot2)
  library(cluster)
  library(splines)
  
  capwords <- function(s, strict = FALSE){
    cap <- function(s){
      paste(toupper(substring(s,1,1)),{s <- substring(s,2); if(strict) tolower(s) else s},
            sep = "", collapse = " ")
      }
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  }
  
  Fun_byIndex<-function(X, index, fun){
    tapply(X, INDEX=index, FUN = fun)
  }
  
  if(is.null(FUNcluster)){
    FUNcluster <- function(x, k, ...){
      clus <- cutree(agnes(x, method=clustering_method, metric=clustering_metric, ...), k=k)
      return(list("cluster"=clus))
    }
  }
  if(!is.function(FUNcluster)){
    stop("the 'FUNcluster' supplied is not a function")
  }
  
  if(is.null(title)){
    if(desc){
      mydesc <- gmt$geneset.descriptions[which(gmt$geneset.names==geneset.name)]
      mytitle <- paste(geneset.name, "\n", mydesc, "\n", indiv, "\n", sep="")
      main.cex <- main.cex*0.15
      lab.cex <- lab.cex*0.5
      axis.cex <- axis.cex*0.5
    }else{
      mytitle <- paste(geneset.name, "\n", indiv, "\n", sep="")
    }
  }else{
    if(title==""){
      mytitle <- NULL
    }else{
      mytitle <- title
    }
  }
  
  interest <- which(gmt$geneset.names==geneset.name)
  if(length(interest)==0){
    stop("The 'geneset.name' supplied is not in the 'gmt'")
  }
  if(is.data.frame(expr)){
    select_probe <- intersect(rownames(expr), unique(gmt$genesets[[interest]]))
    data_sel <- as.matrix(expr[select_probe,])
  }else if(is.list(expr)){
    expr_sel <- expr[[interest]]
    expr_sel <- expr_sel[, , order(as.numeric(dimnames(expr_sel)[[3]]))]
    data_sel <- matrix(expr_sel, nrow=dim(expr_sel)[1], ncol=dim(expr_sel)[2]*dim(expr_sel)[3])
    rownames(data_sel) <- dimnames(expr_sel)[[1]]
    select_probe <- dimnames(expr_sel)[[1]]
    TimePoint <- sort(as.numeric(rep(dimnames(expr_sel)[[3]], dim(expr_sel)[2])))
    Patient_ID <- rep(dimnames(expr_sel)[[2]], dim(expr_sel)[3])
  }
  
  
      
  if(indiv=="genes"){
    data_stand <- t(apply(X=data_sel, MARGIN=1, FUN=scale))
    data_stand_MedianByTP <- t(apply(X=data_stand, MARGIN=1, FUN=Fun_byIndex, index=as.factor(TimePoint), fun=aggreg.fun))
  }else if(indiv=="patients"){
    data_stand <- t(apply(X=data_sel, MARGIN=1, FUN=scale))
    data_tocast<-cbind.data.frame(TimePoint, Patient_ID, "M" = apply(X=data_stand, MARGIN=2, FUN=aggreg.fun))
    data_stand_MedianByTP <- as.matrix(acast(data_tocast, formula="Patient_ID~TimePoint", value.var="M"))
  }
  
  if(!is.null(baseline)){
    colbaseline <- which(sort(unique(TimePoint))==baseline)
    if(length(colbaseline)==0){
      stop("the 'baseline' value used is not one of the time points in 'TimePoint'...\n\n")
    }
    data_stand_MedianByTP <- data_stand_MedianByTP-data_stand_MedianByTP[,colbaseline]
  }

  
  if(clustering | showTrend){
    if(verbose){
      cat("Optimally clustering...\n")
    }
    kmax <- ifelse(dim(data_stand_MedianByTP)[1]>4, max_trends, dim(data_stand_MedianByTP)[1]-1)
    if(kmax>=2){
      cG <- clusGap(x=data_stand_MedianByTP, FUNcluster=FUNcluster, K.max=kmax, B=B, verbose=FALSE)
      nc <- maxSE(f = cG$Tab[, "gap"], SE.f = cG$Tab[, "SE.sim"], method = methodOptiClust)
      clust <- FUNcluster(data_stand_MedianByTP, k=nc)$cluster
    }else{
      nc <- 1
      clust <- rep(1, dim(data_stand_MedianByTP)[1])
    }
    
    medoids <- as.data.frame(t(apply(X=data_stand_MedianByTP, MARGIN=2, FUN=Fun_byIndex, index=clust, fun=trend.fun)))
    if(dim(medoids)[1]==1){
      medoids <- cbind.data.frame("TimePoint"= colnames(medoids), "1"=t(medoids))
    }else{
      medoids <- cbind.data.frame("TimePoint"= rownames(medoids), medoids)
    }
    if(verbose){
        cat("DONE\n")
    }
  }else{
    medoids <- cbind.data.frame("TimePoint"=colnames(data_stand_MedianByTP), "1"='NA')
    clust <- rep(NA, dim(data_stand_MedianByTP)[1])
  }
  meltedData <- melt(cbind.data.frame("Probe_ID"=rownames(data_stand_MedianByTP), "Cluster"=clust, data_stand_MedianByTP), id.vars=c("Probe_ID", "Cluster"), variable.name="TimePoint")
  meltedStats <- melt(medoids, id.vars="TimePoint", variable.name="Cluster")
  meltedData$Cluster <- as.factor(meltedData$Cluster)
  

  meltedData$TimePoint <- paste(time_unit, meltedData$TimePoint, sep="")
  meltedStats$TimePoint <- paste(time_unit, meltedStats$TimePoint, sep="")
  
  if(is.null(y.lim)){
    y.max <- max(abs(meltedData$value))
    y.min <- -y.max
  }else{
    y.max <- y.lim[2]
    y.min <- y.lim[1]
  }
  if(is.null(x.lim)){
    x.lim <- unique(meltedData$TimePoint)
  }
  p <- (ggplot(meltedData, aes(x=TimePoint, y=value)) 
        + geom_hline(aes(y = 0), linetype=1, colour='grey50', size=0.4)
  )
   
  if(!clustering){
     p <- (p
           + geom_line(aes(group=Probe_ID, colour=Probe_ID), size=0.7)
           + scale_colour_manual(guide='none', name='probe ID', values=rainbow(length(select_probe)))
     )
  }else{
    p <- (p
          + geom_line(aes(group=Probe_ID, colour=Cluster), size=0.7)
          + guides(colour = guide_legend(override.aes=list(size=1, fill="white"), keywidth=2*lab.cex, 
                                         title.theme=element_text(size = 15*lab.cex, angle=0),
                                         label.theme=element_text(size = 9*lab.cex, angle=0)
                                         )
                   )
    )
  }

  p <- (p
        + ylab(paste(capwords(aggreg.fun), 'of standardized gene expression'))
        + xlab('Time')
        + ggtitle(mytitle)
        + theme(plot.title=element_text(size = 35*main.cex))
        + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(colour='grey40', fill = 'white'))
        + ylim(y.min, y.max)
        + xlim(x.lim)
        + theme(axis.title.y = element_text(size = 25*lab.cex, angle = y.lab.angle, vjust=0.3), axis.text.y = element_text(size=18*axis.cex, colour = 'grey40')) 
        + theme(axis.title.x = element_text(size = 25*lab.cex, angle = 0, vjust=-0.9), axis.text.x = element_text(size=18*axis.cex, colour = 'grey40', angle=x.axis.angle, vjust=0.5, hjust=0.5))
        + theme(plot.margin=unit(c(0.5, 0.5, 0.7, 1), 'lines'))
        + theme(legend.key=element_rect(fill="white"))
        +gg.add
  )
  
  if(showTrend){
    if(!smooth){
      p <- (p + geom_line(data=meltedStats, aes(x=TimePoint, y=value, group=Cluster, linetype=Cluster), size=4))
    }else{
      p <- (p + stat_smooth(formula=y~poly(x,3), data=meltedStats, aes(x=TimePoint, y=value, group=Cluster, linetype=Cluster), size=4, se=FALSE, method="lm", color="black"))
    }
    p <- p + scale_linetype_manual(name=paste("Cluster", capwords(trend.fun)), values=as.numeric(levels(meltedStats$Cluster))+1, 
                                     guide=guide_legend(override.aes=list(size=1), keywidth=2*lab.cex, 
                                                        title.theme=element_text(size = 15*lab.cex, angle=0),
                                                        label.theme=element_text(size = 9*lab.cex, angle=0)
                                     )
    )
  }
  print(p)
}


  