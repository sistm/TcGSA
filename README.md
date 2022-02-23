
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `TcGSA`

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/TcGSA)](https://cran.r-project.org/package=TcGSA)
[![R-CMD-check](https://github.com/sistm/TcGSA/workflows/R-CMD-check/badge.svg)](https://github.com/sistm/TcGSA/actions)
[![Downloads](https://cranlogs.r-pkg.org/badges/TcGSA?color=blue)](https://www.r-pkg.org/pkg/TcGSA)

## Overview

`TcGSA` is a package which performs *Time-course Gene Set Analysis* from
**microarray data**, and provide nice representations of its results.

On top of the CRAN help pdf-file, the following article explains what
TcGSA is about:

> Hejblum, BP, Skinner, J, & Thiébaut, R (2015). Time-Course Gene Set
> Analysis for Longitudinal Gene Expression Data. *PLOS Computational
> Biology*, **11**(6):e1004310. [\<doi:
> 10.1371/journal.pcbi.1004310>](https://doi.org/10.1371/journal.pcbi.1004310)

## Installation

TcGSA imports the `multtest` package which is not available on
[CRAN](https://cran.r-project.org/), but is available on the
[Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/multtest.html)
repository. Before installing TcGSA, be sure to have this `multtest`
package installed. If not, you can do so by running the following:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("multtest")
```

The easiest way to get `TcGSA` is to install it from
[CRAN](https://cran.r-project.org/package=TcGSA):

``` r
install.packages("TcGSA")
```

or to get the development version from
[GitHub](https://github.com/sistm/TcGSA):

``` r
#install.packages("devtools")
devtools::install_github("sistm/TcGSA")
```

## Microarrays vs RNA-seq

`TcGSA` relies on a Gaussian assumption for the expression data, which
is suitable for normalized microarray data. Due to their count and
heteroskedastic nature, RNA-seq data need to be handled differently and
***TcGSA cannot deal with RNA-seq data***. For RNA-seq data, please have
a look at the [Bioconductor package
`dearseq`](https://bioconductor.org/packages/dearseq/) which
incorporates similar functionality for analyzing RNA-seq data.

– Boris Hejblum
