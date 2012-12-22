summary.TcGSA <-function(object, ...){
  nsignif <- dim(signifLRT.TcGSA(object))[1]
  func_form <- object[["func_form"]]
  separatePatients <- object[["separatePatients"]]
  ntg <- ifelse(is.null(object[["group_var"]]),1,length(levels(object[["group_var"]])))
  ngs <- length(object[["GeneSets_gmt"]]$geneset.name)
  res  <- list("func_form"=func_form, "separatePatients"=separatePatients, "ntg"=ntg, "ngs"=ngs, "nsignif"=nsignif)
  class(res) <- "summary.TcGSA"
  
  return(res)
}