plotPat.TcGSA <-
function(x, threshold=0.05, myproc="BY", nbsimu_pval=1e+06, 
         expr, Patient_ID, TimePoint, 
         baseline=NULL, only.signif=TRUE,
         group.var=NULL, Group_ID_paired=NULL, ref=NULL, group_of_interest=NULL,
         FUNcluster=NULL, clustering_metric="euclidian", clustering_method="ward", B=500,
         max_trends=4, aggreg.fun="median",
         methodOptiClust = "firstSEmax",
         verbose=TRUE,
         clust_trends=NULL,
         N_clusters=NULL, myclusters=NULL, label.clusters=NULL, prev_rowCL=NULL,
         descript=TRUE, plotAll=TRUE,
         color.vec=c("darkred", "#D73027", "#FC8D59", "snow", "#91BFDB", "#4575B4", "darkblue"),
         legend.breaks=NULL,
         label.column=NULL, time_unit="", 
         cex.label.row=1, cex.label.column=1, margins=c(5, 25), 
         heatKey.size=1, dendrogram.size=1, heatmap.height=1, heatmap.width=1,
         cex.clusterKey=1, cex.main=1,
         horiz.clusterKey=TRUE,
         main=NULL, subtitle=NULL, 
         ...){
  
  library(gplots)
  
  Fun_byIndex<-function(X, index, fun){
    tapply(X, INDEX=index, FUN = fun)
  }
  
  gmt <- x[["GeneSets_gmt"]]
  
  if(!is.null(baseline)){
    if(!(baseline %in% unique(TimePoint))){
      stop("The 'baseline' value used is not one of the time points in 'TimePoint'...\n\n")
    }
  }
  
  if(is.null(main)){
    mymain <-paste("Median trends", subtitle, sep="")
  }else if(!is.null(subtitle)){
    mymain <-paste(main, "\n", subtitle, sep="")
  }else{
    mymain <- main
  }
  
  
  
  if(is.null(clust_trends)){
    clust_trends <- clustTrend(x=x, expr=expr, Patient_ID=Patient_ID, TimePoint=TimePoint, baseline=baseline, only.signif=TRUE,
                               group.var=group.var, Group_ID_paired=Group_ID_paired, ref=ref, group_of_interest=group_of_interest,
                               FUNcluster=FUNcluster, clustering_metric=clustering_metric, clustering_method=clustering_method, B=B,
                               max_trends=max_trends, aggreg.fun=aggreg.fun,
                               methodOptiClust = methodOptiClust,
                               indiv="genes",
                               verbose=verbose
    )
  }else if(class(clust_trends)!="ClusteredTrends"){
    stop("The 'clust_trends' argument is not of the class 'ClusteredTrends', see the clustTrend function")
  }
  
  if(verbose){cat("Initializing clustering on all the patients...\n")}
  hc <- plot.TcGSA(x=x, threshold=threshold, myproc=myproc, nbsimu_pval=nbsimu_pval, 
                   expr=expr, Patient_ID=Patient_ID, TimePoint=TimePoint, 
                   baseline=baseline, only.signif=only.signif,
                   group.var=group.var, Group_ID_paired=Group_ID_paired, ref=ref, group_of_interest=group_of_interest,
                   FUNcluster=FUNcluster, clustering_metric=clustering_metric, clustering_method=clustering_method, B=B,
                   max_trends=max_trends, aggreg.fun=aggreg.fun,
                   methodOptiClust = methodOptiClust,
                   indiv="genes",
                   verbose=verbose,
                   clust_trends=clust_trends,
                   N_clusters=N_clusters, myclusters=myclusters, label.clusters=label.clusters, prev_rowCL=prev_rowCL,
                   descript=descript, plot=plotAll,
                   color.vec=color.vec,
                   legend.breaks=legend.breaks,
                   label.column=label.column, time_unit=time_unit, 
                   cex.label.row=cex.label.row, cex.label.column=cex.label.column, margins=margins, heatKey.size=heatKey.size, dendrogram.size=dendrogram.size, heatmap.height=heatmap.height, heatmap.width=heatmap.width,
                   cex.clusterKey=cex.clusterKey, cex.main=cex.main,
                   horiz.clusterKey=horiz.clusterKey,
                   main="Median trends", subtitle="over all patients")

  if(is.null(prev_rowCL)){
    clRows=TRUE
    if(only.signif){
      signif <- multtest.TcGSA(x, threshold, myproc, nbsimu_pval)
      select <- which(signif$adj_pval<0.05)
      if(!length(select)>0){
        stop("No gene sets significant")
      }
      subtitle <- paste(subtitle, "\n", length(which(signif$adj_pval<0.05)), "/", length(signif$adj_pval), " gene sets significant with ", x[["func_form"]], " shape", sep="")
    }else{
      select <- 1:length(gmt$geneset.names)
      paste(subtitle, "\n", length(signif$adj_pval), " gene sets", sep="")
    }
  }else{
    if(!is.null(prev_rowCL$ddr)){
      clRows=prev_rowCL$ddr
    }else{
      clRows=as.dendrogram(prev_rowCL)
    }
    
    select <- match(prev_rowCL$labels, gmt$geneset.names)
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
  
  pat <- levels(Patient_ID)
  #par('mfrow'=c(ceiling(length(pat)/4),4))
  

  for (i in 1:length(pat)){
    p <- pat[i]
    clust_trends_p <- clust_trends
    
    if(verbose){
      cat(paste("Patient ", i, "/", length(pat),":", sep=""))
    }
    
    for(interest in select){
      gs <- gmt$geneset.names[interest]
      if(is.null(group.var)){
        if(is.data.frame(expr)){
          select_probe <- intersect(rownames(expr), unique(gmt$genesets[[interest]]))
          data_sel <- as.matrix(expr[select_probe, ])
        }else if(is.list(expr)){
          expr_sel <- expr[[interest]]
          expr_sel <- expr_sel[, , order(as.numeric(dimnames(expr_sel)[[3]]))]
          data_sel <- matrix(expr_sel, nrow=dim(expr_sel)[1], ncol=dim(expr_sel)[2]*dim(expr_sel)[3])
          select_probe <- dimnames(expr_sel)[[1]]
          rownames(data_sel) <- select_probe
          TimePoint <- sort(as.numeric(rep(dimnames(expr_sel)[[3]], dim(expr_sel)[2])))
          Patient_ID <- rep(dimnames(expr_sel)[[2]], dim(expr_sel)[3])
        }
        
        data_stand <- t(apply(X=data_sel, MARGIN=1, FUN=scale))
        # So the genes expression is comparable over the patients.
        
      }else{
        if(!is.null(baseline)){
          stop("the 'baseline' argument is not NULL while a grouping variable is supplied in 'group.var'...\n")
        }
        if(is.data.frame(expr)){
          select_probe <- intersect(rownames(expr), unique(gmt$genesets[[interest]]))
          data_sel <- as.matrix(expr[select_probe,])
        }else if(is.list(expr)){
          expr_sel <- expr[[interest]]
          expr_sel <- expr_sel[, , order(as.numeric(dimnames(expr_sel)[[3]]))]
          data_sel <- matrix(expr_sel, nrow=dim(expr_sel)[1], ncol=dim(expr_sel)[2]*dim(expr_sel)[3])
          select_probe <- dimnames(expr_sel)[[1]]
          rownames(data_sel) <- select_probe
          if(!is.null(Group_ID_paired)){
            Group_ID_paired <- Group_ID_paired[order(TimePoint)] # watch out for the ordering
          }
          TimePoint <- sort(as.numeric(rep(dimnames(expr_sel)[[3]], dim(expr_sel)[2])))
          Patient_ID <- rep(dimnames(expr_sel)[[2]], dim(expr_sel)[3])
        }
        data_stand_ref <- t(apply(X=data_sel[,group.var==ref], MARGIN=1, FUN=scale))
        data_stand_interest <- t(apply(X=data_sel[,group.var==group_of_interest], MARGIN=1, FUN=scale))
        
        if(is.null(Group_ID_paired)){
          data_stand_ByTP_ref <- t(apply(X=data_stand_ref, MARGIN=1, FUN=Fun_byIndex, index=as.factor(TimePoint)[group.var==ref], fun=aggreg.fun))
          data_stand_ByTP_interest <- t(apply(X=data_stand_interest, MARGIN=1, FUN=Fun_byIndex, index=as.factor(TimePoint)[group.var==group_of_interest], fun=aggreg.fun))                         
          data_stand_ByTP <- data_stand_ByTP_interest-data_stand_ByTP_ref
        }else{
          data_diff <- t(apply(X=cbind.data.frame(data_stand_interest, -data_stand_ref), MARGIN=1, FUN=Fun_byIndex, 
                                index=(as.factor(c(TimePoint[group.var==group_of_interest], TimePoint[group.var==ref])):as.factor(c(as.character(Group_ID_paired)[group.var==group_of_interest], as.character(Group_ID_paired)[group.var==ref]))),
                                fun=sum))
          data_stand <- t(apply(X=data_diff, MARGIN=1, FUN=Fun_byIndex, index=sort(as.factor(TimePoint[group.var==group_of_interest])), fun=aggreg.fun))
        }
      }
      
      
      data_stand_ByTP <- data_stand[,Patient_ID==p]
      
      if(!is.null(baseline)){
        colbaseline <- which(sort(unique(TimePoint))==baseline)
        if(length(colbaseline)==0){
          stop("the 'baseline' value used is not one of the time points in 'TimePoint'...\n")
        }
        data_stand_ByTP <- data_stand_ByTP-data_stand_ByTP[,colbaseline]
      }
      
      medoids <- as.data.frame(t(apply(X=data_stand_ByTP, MARGIN=2, FUN=Fun_byIndex, index=clust_trends_p$GenesPartition[[gs]], fun="median")))
      if(dim(medoids)[1]==1){
        medoids <- cbind.data.frame("TimePoint"= TimePoint[Patient_ID==p], "1"=t(medoids))
        rownames(medoids) <- TimePoint[Patient_ID==p]
      }else{
        medoids <- cbind.data.frame("TimePoint"= rownames(medoids), medoids)
      }
      clust_trends_p$ClustMeds[[gs]] <- medoids
    }

    plot.TcGSA(x=x, 
               threshold=threshold, myproc=myproc, nbsimu_pval=nbsimu_pval, 
               expr=NULL, Patient_ID=NULL, TimePoint=TimePoint, 
               baseline=baseline, only.signif=only.signif,
               group.var=group.var, Group_ID_paired=Group_ID_paired, ref=ref, group_of_interest=group_of_interest,
               FUNcluster=FUNcluster, clustering_metric=clustering_metric, clustering_method=clustering_method, B=B,
               max_trends=max_trends, aggreg.fun=aggreg.fun,
               methodOptiClust = methodOptiClust,
               indiv="genes",
               verbose=verbose,
               clust_trends=clust_trends_p,
               N_clusters=NULL, myclusters=hc$myclusters, label.clusters=label.clusters, prev_rowCL=hc,
               descript=descript, plot=TRUE,
               color.vec=color.vec,
               legend.breaks=hc$legend.breaks,
               label.column=label.column, time_unit=time_unit, 
               cex.label.row=cex.label.row, cex.label.column=cex.label.column, margins=margins, heatKey.size=heatKey.size, dendrogram.size=dendrogram.size, heatmap.height=heatmap.height, heatmap.width=heatmap.width,
               cex.clusterKey=cex.clusterKey, cex.main=cex.main,
               horiz.clusterKey=horiz.clusterKey,
               main="Median Trends", subtitle=paste("Patient", p)
    )
      
    if(verbose){
      cat(" done\n")
    }           
  }
  return(hc)
}
