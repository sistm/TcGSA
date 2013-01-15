clustTrend <- 
function(x,
         expr, Patient_ID, TimePoint, baseline=NULL, only.signif=TRUE,
         group.var=NULL, Group_ID_paired=NULL, ref=NULL, group_of_interest=NULL,
         FUNcluster=NULL, clustering_metric="euclidian", clustering_method="ward", B=100,
         max_trends=4, aggreg.fun="median", trend.fun="median",
         methodOptiClust = "firstSEmax",
         indiv="genes",
         verbose=TRUE
         ){
#  library(cluster)
  Fun_byIndex<-function(X, index, fun, ...){
    tapply(X, INDEX=index, FUN = fun, ...)
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
  
  gmt <- x[["GeneSets_gmt"]]
  if(only.signif){
    GSsig <- signifLRT.TcGSA(x)
    GeneSetsList <- GSsig$GeneSet
  }
  else{
    GeneSetsList <- gmt$geneset.names
  }
  
  if(!is.null(group.var)){
    if(is.null(ref)){
      ref <- levels(group.var)[1]
    }
    if(is.null(group_of_interest)){
      group_of_interest <- levels(group.var)[2]
    }
  }
  
  NbClust <- numeric(length(GeneSetsList))
  names(NbClust) <- GeneSetsList
  ClustsMeds <- list()
  GenesPartition <- list()
  for(gs in GeneSetsList){
    interest <- which(gmt$geneset.names==as.character(gs))
    if(is.null(group.var)){
      if(is.data.frame(expr)){
        select_probe <- intersect(rownames(expr), unique(gmt$genesets[[interest]]))
        data_sel <- as.matrix(expr[select_probe,])
      }
      else if(is.list(expr)){
        expr_sel <- expr[[interest]]
        expr_sel <- expr_sel[, , order(as.numeric(dimnames(expr_sel)[[3]]))]
        # data_sel <- matrix(expr_sel, nrow=dim(expr_sel)[1], ncol=dim(expr_sel)[2]*dim(expr_sel)[3])
        data_sel <- acast(melt(expr_sel, varnames=c("Probe_ID", "Patient_ID", "TimePoint")), 
        			formula=Probe_ID ~ TimePoint + Patient_ID)
        # rownames(data_sel) <- dimnames(expr_sel)[[1]]
        select_probe <- dimnames(expr_sel)[[1]]
        TimePoint <- sort(as.numeric(rep(dimnames(expr_sel)[[3]], dim(expr_sel)[2])))
        Patient_ID <- rep(dimnames(expr_sel)[[2]], dim(expr_sel)[3])
      }
      
      data_stand <- t(apply(X=data_sel, MARGIN=1, FUN=scale))
      if(indiv=="genes"){
        data_stand_ByTP <- t(apply(X=data_stand, MARGIN=1, FUN=Fun_byIndex, index=as.factor(TimePoint), fun=aggreg.fun, na.rm=T))
      }
      else if(indiv=="patients"){
        data_tocast<-cbind.data.frame(TimePoint, Patient_ID, "M" = apply(X=data_stand, MARGIN=2, FUN=aggreg.fun, na.rm=T))
        data_stand_ByTP <- as.matrix(acast(data_tocast, formula="Patient_ID~TimePoint", value.var="M"))
      }
      
      if(!is.null(baseline)){
        colbaseline <- which(sort(unique(TimePoint))==baseline)
        if(length(colbaseline)==0){
          stop("the 'baseline' value used is not one of the time points in 'TimePoint'...\n")
        }
        data_stand_ByTP <- data_stand_ByTP-data_stand_ByTP[,colbaseline]
      }
      
    }
    else{
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
        rownames(data_sel) <- dimnames(expr_sel)[[1]]
        select_probe <- dimnames(expr_sel)[[1]]
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
        data_stand <- t(apply(X=cbind.data.frame(data_stand_interest, -data_stand_ref), MARGIN=1, FUN=Fun_byIndex, 
                            index=(as.factor(c(TimePoint[group.var==group_of_interest], TimePoint[group.var==ref])):as.factor(c(as.character(Group_ID_paired)[group.var==group_of_interest], as.character(Group_ID_paired)[group.var==ref]))),
                            fun=sum))
        data_stand_ByTP <- t(apply(X=data_stand, MARGIN=1, FUN=Fun_byIndex, index=sort(as.factor(TimePoint[group.var==group_of_interest])), fun=aggreg.fun))
      }
    }
   
    kmax <- ifelse(dim(data_stand_ByTP)[1]>4, max_trends, dim(data_stand_ByTP)[1]-1)
    if(kmax>=2){
      cG <- clusGap(x=data_stand_ByTP, FUNcluster=FUNcluster, K.max=kmax, B=B, verbose=FALSE)
      nc <- maxSE(f = cG$Tab[, "gap"], SE.f = cG$Tab[, "SE.sim"], method = methodOptiClust)
      clust <- FUNcluster(data_stand_ByTP, k=nc)$cluster
    }else{
      nc <- 1
      clust <- rep(1, dim(data_stand_ByTP)[1])
    }
    
    medoids <- as.data.frame(t(apply(X=data_stand_ByTP, MARGIN=2, FUN=Fun_byIndex, index=clust, fun=trend.fun)))
    if(dim(medoids)[1]==1){
      medoids <- cbind.data.frame("TimePoint"= colnames(medoids), "1"=t(medoids))
    }else{
      medoids <- cbind.data.frame("TimePoint"= rownames(medoids), medoids)
    }
    
    NbClust[gs] <- nc
    ClustsMeds[[gs]] <- medoids
    GenesPartition[[gs]] <- clust
    names(GenesPartition[[gs]]) <- select_probe
    if(verbose){
      cat(paste(which(GeneSetsList==gs), "/", length(GeneSetsList), " gene sets clustered\n", sep=""))
    }
  }
  res <- list("NbClust"=NbClust, "ClustMeds"=ClustsMeds, "GenesPartition"=GenesPartition, "MaxNbClust"=max_trends)
  class(res) <- "ClusteredTrends"
  return(res)
}
