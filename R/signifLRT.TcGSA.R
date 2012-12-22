signifLRT.TcGSA <-
function(tcgsa, threshold=0.05, myproc="BY", nbsimu_pval = 1e+06, write=F, txtfilename=NULL, directory=NULL){  
  gmt  <-  tcgsa[["GeneSets_gmt"]]
  signif <- multtest.TcGSA(tcgsa, threshold, myproc, nbsimu_pval)

  signif_mod <- gmt$geneset.name[which(signif$adj_pval<threshold)]
  if(!is.null(signif_mod)){
    signif_desc <- gmt$geneset.descriptions[which(signif$adj_pval<threshold)]
    AdjPval <- signif$adj_pval[which(signif$adj_pval<threshold)]
  }else{
    signif_desc <- NULL
    AdjPval <- NULL
  }

  Res_Linear_Mod_FDR <- cbind.data.frame("GeneSet"=signif_mod, "AdjPval"=AdjPval, "desc"=signif_desc)
  if(dim(Res_Linear_Mod_FDR)[1]>0){
    Res_Linear_Mod_FDR <- Res_Linear_Mod_FDR[mixedorder(gsub('.',"a",as.character(Res_Linear_Mod_FDR$GeneSet), fixed=TRUE)),]
  }
  
 if(write){
   if(!is.null(txtfilename)){
     if (is.null(directory)){
        directory <- getwd()
        cat("Warning: 'directory' argument is empty, output file written in the current working directory")
      }
      write.table(Res_Linear_Mod_FDR, file=paste(directory, txtfilename, sep="/"), row.names=FALSE, sep="\t")
    }else{
      cat("ERROR: could not write the significant results file because the argument 'txtfilename' is empty")
    }
  }
  
  return(Res_Linear_Mod_FDR)
}
