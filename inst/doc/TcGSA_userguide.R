## ----pre, echo=FALSE, warning=FALSE, include=FALSE-----------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE, eval=TRUE, cache=TRUE)
#rmarkdown::render("vignettes/TcGSA_userguide.Rmd")

## ----GEOquery, include=FALSE, message=FALSE, cache=FALSE-----------------
if (!requireNamespace("GEOquery", quietly = TRUE)) {
	source("https://bioconductor.org/biocLite.R")
	biocLite("GEOquery")
}

