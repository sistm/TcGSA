\name{estim.TcGSA}
\alias{estim.TcGSA}

\title{
Computing Genes Dynamics out of Estimations of Effects in Gene Sets Modelling
}

\description{
This function computes genes dynamics inside gene sets out of the estimations obtained from mixed models in \code{\link{TcGSA.LR}}.
}
\usage{
estim.TcGSA(tcgsa, expr, Patient_ID, TimePoint)
}

\arguments{
  \item{tcgsa}{
a TcGSA object
  }
    \item{expr}{
a matrix or dataframe of gene expression.  Its dimension are \eqn{n}x\eqn{p}, with the \eqn{p} samples in column and the \eqn{n} genes in row.
  }
  \item{Patient_ID}{
a factor of length \eqn{p} that is in the same order as the columns of \code{expr} and that contains the patient identifier of each sample.%  See Details.
  }
  \item{TimePoint}{
a numeric vector or a factor of length \eqn{p} that is in the same order as \code{TimePoint}, and that contains the time points at which gene expression was measured.%  See Details.
  }
}

%\details{
%TODO
%}

\value{
\code{estim.TcGSA} returns a list of 3 dimensions arrays.  Each element of the list (i.e. each array) corresponds to the estimations of gene expression dynamics for each of the gene sets under scrutiny (obtained from mixed models).  The first dimension of those arrays is the genes included in the concerned gene set, the second dimension is the \code{Patient_ID}, and the third dimension is the \code{TimePoint}.  The values inside those arrays are estimated gene expressions.
}

\references{
Hejblum, B.P., Skinner, J., Thiebaut, R., 2012, TcGSA: Time-course gene set analysis of longitudinal gene expression, with an application to a therapeutic HIV vaccine trial, \bold{submitted}.
}

\author{
Boris P. Hejblum
}

\seealso{
  \code{\link{TcGSA.LR}}, \code{\link{TcGSA.LR.parallel}}
}
\examples{
data(data_simu_TcGSA)

tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, Patient_ID=design$Patient_ID,
                           TimePoint=design$TimePoint, func = "linear")
genesDyn_1grp <- estim.TcGSA(tcgsa=tcgsa_sim_1grp, expr=expr_1grp, Patient_ID=design$Patient_ID, 
                             TimePoint=design$TimePoint)

summary(tcgsa_sim_1grp)
geneset=1
timepoint="1"
genesDyn_1grp[[geneset]][, ,timepoint]

geneset=4
gene="167"
genesDyn_1grp[[geneset]][gene, ,]

geneset=5
patient="P9"
genesDyn_1grp[[geneset]][,patient ,]

plot1GS(expr=genesDyn_1grp, TimePoint=design$TimePoint, Patient_ID=design$Patient_ID, gmt=gmt_sim,
        geneset.name="Gene set 3",
        indiv="genes",
        time_unit="H")

}
