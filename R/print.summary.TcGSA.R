print.summary.TcGSA <-function(x, ...){
  cat("\t\tA TcGSA object")
  cat("\n")
  cat("Form of the time trend:")
  cat("\n\t")
  cat(x[["func_form"]])
  cat("\n")
  cat("Number of treatment groups:")
  cat("\n\t")
  cat(x[["ntg"]])
  cat("\n")
  if(x[["separatePatients"]]){
    cat("Number of gene sets tested for discriminating time trends among patients:") 
  }else{
    cat("Number of gene sets tested for significant time trend:")
  }
  cat("\n\t")
  cat(x[["ngs"]])
  cat("\n\n")
  cat("Number of significant gene sets at a 5% FDR (Benjamini & Yekutieli step-up procedure):")
  cat("\n\t")
  cat(x[["nsignif"]])
  cat(" out of ")
  cat(x[["ngs"]])
  cat(" gene sets")
}