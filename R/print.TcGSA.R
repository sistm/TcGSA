print.TcGSA <- function(x, ...){
  cat("\t\tA TcGSA object")
  cat("\n")
  cat("Form of the time trend:")
  cat("\n\t")
  cat(x[["func_form"]])
  cat("\n")
  cat("Number of treatment groups:")
  cat("\n\t")
  cat(ifelse(is.null(x[["group_var"]]),1,length(levels(x[["group_var"]]))))
  cat("\n")
  if(x[["separatePatients"]]){
    cat("Number of gene sets tested for discriminating time trends among patients:") 
  }else{
    cat("Number of gene sets tested for significant time trend:")
  }
  cat("\n\t")
  cat(length(x[["GeneSets_gmt"]]$geneset.name))
}