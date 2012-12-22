\name{signifLRT.TcGSA}
\alias{signifLRT.TcGSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifiying the Significant Gene Sets
}
\description{
A function that identifies the significant gene sets in an object of class '\code{TcGSA}'.
}
\usage{
signifLRT.TcGSA(tcgsa, threshold = 0.05, myproc = "BY", nbsimu_pval = 1e+06,
                write = F, txtfilename = NULL, directory = NULL)
}

\arguments{
  \item{tcgsa}{
a \code{tcgsa} object.
}
  \item{threshold}{
the threshold at which the FDR or the FWER should be controlled.
}
  \item{myproc}{
  a vector of character strings containing the names of the multiple testing procedures for which adjusted p-values are to be computed.  This vector should include any of the following: "\code{Bonferroni}", "\code{Holm}", "\code{Hochberg}", "\code{SidakSS}", "\code{SidakSD}", "\code{BH}", "\code{BY}", "\code{ABH}", "\code{TSBH}". See \code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}} for details.  Default is "\code{BY}", the Benjamini & Yekutieli (2001) step-up FDR-controlling procedure (general dependency structures).  In order to control the FWER(in case of an analysis that is more a hypothesis confirmation than an exploration of the expression data), we recommand to use "\code{Holm}", the Holm (1979) step-down adjusted p-values for strong control of the FWER.
}
  \item{nbsimu_pval}{
the number of observations under the null distribution to be generated in order to compute the p-values. Default is \code{1e+06}
}
  \item{write}{
logical flag enabling the export of the results as a table in a .txt file.  Default is \code{FALSE}.
}
  \item{txtfilename}{
a character string with the name of the .txt file in which the results table is to be written, if \code{write} is \code{TRUE}. Default is \code{NULL}.
}
  \item{directory}{
if \code{write} is \code{TRUE}, a character string with the directory of the .txt file in which the results table is to be written, if \code{write} is \code{TRUE}. Default is \code{NULL}.
}
}

\value{
\code{signifLRT.TcGSA} returns a data frame with \eqn{p} rows (one row for each significant gene set) and the 3 following variables:
\item{GeneSet}{the significant gene set name from the gmt object.}
\item{AdjPval}{the adjusted p-value corresponding to the signicant gene set.}
\item{desc}{the significant gene set description from the gmt object.}
}

\references{
Hejblum, B.P., Skinner, J., Thiebaut, R., 2012, TcGSA: Time-course gene set analysis of longitudinal gene expression, with an application to a therapeutic HIV vaccine trial, \bold{submitted}.
}

\author{
Boris P. Hejblum
}

\seealso{
\code{\link{multtest.TcGSA}}, \code{\link{TcGSA.LR}}
}
\examples{
data(data_simu_TcGSA)
tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, Patient_ID=design$Patient_ID,
                           TimePoint=design$TimePoint, func = "linear")
sgnifs <- signifLRT.TcGSA(tcgsa_sim_1grp, threshold = 0.05, myproc = "BY",
                          nbsimu_pval = 1000, write=FALSE)
sgnifs
}
