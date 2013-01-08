plot.TcGSA <-
  function(x, threshold=0.05, myproc="BY", nbsimu_pval=1e+06, 
           expr, Patient_ID, TimePoint, 
           baseline=NULL, only.signif=TRUE,
           group.var=NULL, Group_ID_paired=NULL, ref=NULL, group_of_interest=NULL,
           FUNcluster=NULL, clustering_metric="euclidian", clustering_method="ward", B=500,
           max_trends=4, aggreg.fun="median",
           methodOptiClust = "firstSEmax",
           indiv="genes",
           verbose=TRUE,
           clust_trends=NULL,
           N_clusters=NULL, myclusters=NULL, label.clusters=NULL, prev_rowCL=NULL,
           descript=TRUE, plot=TRUE,
           color.vec=c("darkred", "#D73027", "#FC8D59", "snow", "#91BFDB", "#4575B4", "darkblue"),
           legend.breaks=NULL,
           label.column=NULL, time_unit="", 
           cex.label.row=1, cex.label.column=1, margins=c(5, 25), heatKey.size=1, dendrogram.size=1, heatmap.height=1, heatmap.width=1,
           cex.clusterKey=1, cex.main=1,
           horiz.clusterKey=TRUE,
           main=NULL, subtitle=NULL, 
           ...){
#		library(gplots)
    gmt <- x[["GeneSets_gmt"]]
    
    if(!is.null(baseline)){
      if(!(baseline %in% unique(TimePoint))){
        stop("The 'baseline' value used is not one of the time points in 'TimePoint'...\n\n")
      }
    }
    
    if(is.null(main)){
      if(!is.null(subtitle)){
        mymain <-paste("Median trends", "\n", subtitle, sep="")
      }else{
        mymain <-paste("Median trends\n", "over all patients", sep="")
      }
    }else {
      if(!is.null(subtitle)){
        mymain <-paste(main, "\n", subtitle, sep="")
      }else{
        mymain <- main
      }
    }
    
    if(is.null(prev_rowCL)){
      clRows=TRUE
      signif <- multtest.TcGSA(x, threshold, myproc, nbsimu_pval)
      select <- which(signif$adj_pval<0.05)
      subtitle <- paste(subtitle, "\n", length(which(signif$adj_pval<0.05)), "/", length(signif$adj_pval), " gene sets significant with ", x[["func_form"]], " shape", sep="")
    }else{
      if(!is.null(prev_rowCL$ddr)){
        clRows=prev_rowCL$ddr
      }else{
        clRows=as.dendrogram(prev_rowCL)
      }
      
      select <- match(prev_rowCL$geneset.names, gmt$geneset.names)
      if(length(which(is.na(select)))>0){
        select <- match(prev_rowCL$labels, gsub(": Undetermined","", paste(gmt$geneset.names, ": ", gmt$geneset.description, sep="")))
      }
      if(length(which(is.na(select)))>0){
        stop("Geneset names used in the previous clustering don't match the 'geneset.names' from the 'gmt' element of the 'x' argument")
      }
      
      if(is.null(myclusters) && !is.null(prev_rowCL$myclusters)){
        myclusters <- prev_rowCL$myclusters
      }  
    }
    
    if(only.signif & is.null(clustTrend)){
      if(!length(select)>0){
        stop("No gene sets significant")
      }
      gmt <- list("genesets"=gmt$genesets[select], "geneset.names"=gmt$geneset.names[select], "geneset.descriptions"=gmt$geneset.descriptions[select])
      class(gmt) = "GSA.genesets"
    }
    
    if(is.null(clust_trends)){
      clust_trends <- clustTrend(x=x, expr=expr, Patient_ID=Patient_ID, TimePoint=TimePoint, baseline=baseline, only.signif=TRUE,
                                 group.var=group.var, Group_ID_paired=Group_ID_paired, ref=ref, group_of_interest=group_of_interest,
                                 FUNcluster=FUNcluster, clustering_metric=clustering_metric, clustering_method=clustering_method, B=B,
                                 max_trends=max_trends, aggreg.fun=aggreg.fun,
                                 methodOptiClust = methodOptiClust,
                                 indiv=indiv,
                                 verbose=verbose
                                 )
    }else if(class(clust_trends)!="ClusteredTrends"){
      stop("The 'clust_trends' argument is not of the class 'ClusteredTrends', see the clustTrend function")
    }
    
    medoids2clust <- acast(melt(clust_trends[["ClustMeds"]], variable.name="Cluster", id.vars="TimePoint"),
                   formula="L1 + Cluster~ TimePoint", value.var="value")
    gsNames <- gsub("_.*$", "", rownames(medoids2clust))
    ncl <- gsub("^.*?_", "", rownames(medoids2clust))
    medoids2clust <- medoids2clust[,order(as.numeric(colnames(medoids2clust)))]
    
    if(!descript){
      rownames(medoids2clust) <- paste(gsub("_", " ", rownames(medoids2clust)), clust_trends[[1]][gsNames], sep="/")
    }else{
      rownames(medoids2clust) <- paste(gsub(": Undetermined","", paste(gmt$geneset.names[match(gsNames, gmt$geneset.names)], ": ", gmt$geneset.description[match(gsNames, gmt$geneset.names)], sep="")),
                                       " ", ncl, "/", clust_trends[[1]][gsNames],
                                       sep="")
    }
    
    map2heat <- medoids2clust
    
    
    if(is.null(prev_rowCL)){
      hc <- hclust(d=dist(map2heat, method = "euclidean"), method="ward")
      row_wt <- rowMeans(x=map2heat, na.rm = TRUE)
      ddr <- reorder(x=as.dendrogram(hc), wts=row_wt)
    }else{
      ddr <- prev_rowCL$ddr
      hc <- prev_rowCL
    }
    
    
    
    if(!is.null(N_clusters) && is.null(myclusters)){
      myclusters <- as.factor(cutree(hc, k=N_clusters))
      myclusters_num <- levels(myclusters)
      if(N_clusters<9){
        levels(myclusters) <- (c(hsv(0.56, 0.9, 1), hsv(0, 0.27, 1), hsv(0.52, 1, 0.5), 
                                hsv(0.12, 0.55, 0.97), hsv(0.83, 0.81, 0.55), hsv(0.66, 0.15, 1),
                                hsv(0.7, 1, 0.7), hsv(0.42, 0.33, 1)
                               )[1:N_clusters]
                              )
      }else{
        levels(myclusters) <- rainbow(N_clusters, start=0.1, end=0.9)
      }
      myclusters <- as.character(myclusters)
    }else if(!is.null(myclusters)){
      myclusters_num <- 1:length(levels(as.factor(myclusters)))
      N_clusters <- length(levels(as.factor(myclusters)))
    }
    
    if(plot){
      
      myhclustward<- function(d, method = "ward", members=NULL){
        hclust(d, method = "ward", members=NULL)
      }
      
      if(is.null(N_clusters)){
        heatKey.size <- 2.6*heatKey.size
        dendrogram.size <- 2*dendrogram.size
      }else{
        heatKey.size <- 0.9*heatKey.size
        dendrogram.size <- 0.8*dendrogram.size
      }
      
      #d <- length(unique(map2heat))
      #if(floor(d/2)!=d/2){d <- d-1}
      maxAbs <- max(abs(min(map2heat)), abs(max(map2heat)))
      if(is.null(legend.breaks)){
        legend.breaks <- c(seq(from=-ceiling(maxAbs*10)/10, to=0, by=ceiling(maxAbs)/100),
                           seq(from=0, to=ceiling(maxAbs*10)/10, by=ceiling(maxAbs)/100)[-1]
        )
        #legend.breaks <- c(seq(from=-maxAbs, to=0, length.out=d/2+1),
        #                   seq(from=0, to=maxAbs, length.out=d/2+1)[-1]
        #)
      }
      
      colnames(map2heat) <- paste(time_unit, colnames(map2heat), sep="")
      try(
        if(!is.null(myclusters)){
          MYheatmap.2(x = map2heat, 
                      Rowv=clRows,
                      Colv=FALSE,
                      hclustfun=myhclustward,
                      dendrogram='row',
                      scale="none",
                      col=colorRampPalette(rev(color.vec))(length(legend.breaks)-1),
                      breaks=legend.breaks,
                      symkey=TRUE,
                      colsep=NULL,
                      rowsep=NULL, #c(1:dim(map2heat)[1]),
                      sepwidth=c(0.01,0.01),
                      trace="none",
                      RowSideColors = myclusters,
                      density.info="none",
                      lmat=matrix(c(5,6,4,3,1,2), nrow=2,ncol=3,byrow=TRUE),
                      lhei=c(0.115*heatKey.size, 0.3*heatmap.height),
                      lwid=c(0.1*dendrogram.size,0.01,0.4*heatmap.width),
                      cexRow = 0.1*cex.label.row + 0.5*cex.label.row/log10(dim(map2heat)[1]),
                      cexCol = 0.1*cex.label.column + 0.5*cex.label.column/log10(dim(map2heat)[2]),
                      margins=margins,
                      main=list(mymain, cex=1.3*cex.main),
                      ...
          )
          if(N_clusters<9){
            legFill<- (c(hsv(0.56, 0.9, 1), hsv(0, 0.27, 1), hsv(0.52, 1, 0.5), 
                         hsv(0.12, 0.55, 0.97), hsv(0.83, 0.81, 0.55), hsv(0.66, 0.15, 1),
                         hsv(0.7, 1, 0.7), hsv(0.42, 0.33, 1)
                        )[1:N_clusters]
                       )
          }else{
            legFill <- rainbow(N_clusters, start=0.1, end=0.9)
            legFill <- legFill[order(legFill)]
          }
          if(is.null(label.clusters) | length(label.clusters)!=N_clusters){
            legend("topright",legend=myclusters_num, fill=legFill, cex=0.7*cex.clusterKey, horiz=horiz.clusterKey)
          }else{
            legend("topright",legend=label.clusters, fill=legFill, cex=0.7*cex.clusterKey, horiz=horiz.clusterKey)
          }        
        }else{
          MYheatmap.2(x = map2heat,
                      Rowv=clRows,
                      Colv=FALSE,
                      hclustfun=myhclustward,
                      dendrogram='row',
                      scale="none",
                      col=colorRampPalette(rev(color.vec))(length(legend.breaks)-1),
                      breaks=legend.breaks,
                      symkey=TRUE,
                      colsep=NULL,
                      rowsep=NULL, #c(1:dim(map2heat)[1]),
                      sepwidth=c(0.01,0.01),
                      trace="none",
                      density.info="none",
                      lmat=matrix(c(4,3,2,1), nrow=2,ncol=2,byrow=TRUE),
                      lhei=c(0.115*heatKey.size,1),
                      lwid=c(0.1*dendrogram.size,1),
                      cexRow = 0.1*cex.label.row + 0.5*cex.label.row/log10(dim(map2heat)[1]),
                      cexCol = 0.1*cex.label.column + 0.5*cex.label.column/log10(dim(map2heat)[2]),
                      margins=margins,
                      main=list(mymain, cex=1.3*cex.main),
                      ...
          )
        }
      )
    }
    hc$legend.breaks <- legend.breaks
    hc$myclusters <- myclusters
    hc$ddr <- ddr
    hc$geneset.names <- gmt$geneset.names[select]
    hc$clust.trends <- clust_trends
    
    if(!is.null(myclusters)){
      clusters <- as.factor(myclusters)
      if(!is.null(label.clusters)){
        levels(clusters) <- label.clusters
      }else{
        levels(clusters) <- myclusters_num
      }
      GSs <- rownames(map2heat)
      clustersExport <- cbind.data.frame("GeneSetTrend"=GSs, "Cluster"=as.character(clusters))
      rownames(clustersExport) <- GSs
      clustersExport  <-  clustersExport[order(clusters),]
    }else{
      clustersExport<- NULL
    }
    
    hc$clusterExport <- clustersExport
    
    return(hc)
  }
