multtest.TcGSA <-
function(tcgsa, threshold=0.05, myproc="BY", nbsimu_pval = 1000000){
  emp <- tcgsa[["fit"]]
  func <- tcgsa[["func_form"]]
  group.var <- tcgsa[["group.var"]]
  separatePatients <- tcgsa[["separatePatients"]]
	splines_DF <- tcgsa[["splines_DF"]]
  
  if(is.null(group.var)){
    if(!separatePatients){
      if(func=="linear"){
        theodist <- c(rchisq(nbsimu_pval/2,df=1), rchisq(nbsimu_pval/2,df=2) )
      }else if(func=="cubic"){
        theodist <- rmixchisq(nbsimu_pval,3,3)
      }else if(func=="splines"){
        theodist <-rmixchisq(nbsimu_pval,splines_DF,splines_DF)
      }
    }else{
      if(func=="linear"){
        theodist <- c(rchisq(nbsimu_pval/2,df=0), rchisq(nbsimu_pval/2,df=1) )
      }else if(func=="cubic"){
        theodist <- rmixchisq(nbsimu_pval,0,3)
      }else if(func=="splines"){
        theodist <-rmixchisq(nbsimu_pval,0,splines_DF)
      }
    }
  }else{
    nbgp <- length(levels(group.var))
    if(func=="linear"){
      theodist <- rmixchisq(nbsimu_pval,((1*nbgp-1)+1),(1*nbgp))
    }else if(func=="cubic"){
      theodist <- rmixchisq(nbsimu_pval,((3*nbgp-1)+1),(3*nbgp))
    }else if(func=="splines"){
      theodist <-rmixchisq(nbsimu_pval,((splines_DF*nbgp-1)+1),splines_DF*nbgp)
    }
  }
 
  emp$raw_pval <- unlist(lapply(emp$LR, FUN=pval_simu, theo_dist=theodist))
  adj_pval <- mt.rawp2adjp(emp$raw_pval, proc=c(myproc),alpha=threshold)
  emp$adj_pval <- adj_pval$adjp[order(adj_pval$index),2]

  return(emp)
}
