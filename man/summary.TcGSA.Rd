\name{summary.TcGSA}
\alias{summary.TcGSA}
\alias{print.summary.TcGSA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing TcGSA
}
\description{
\code{summary} method for class '\code{TcGSA}'
}
\usage{
\S3method{summary}{TcGSA}(object, ...)

\S3method{print}{summary.TcGSA}(x, ...)
}

\arguments{
  \item{object}{
    an object of class'\code{TcGSA}'.
  }  
  \item{x}{
    an object of class'\code{summary.TcGSA}'.
  }
  \item{\dots}{
    further arguments passed to or from other methods.
  }
}

\value{
The function summary.TcGSA returns a list with the following components (list elements):
\item{func_form}{the chosen form for the time trend.}
\item{separatePatients}{a logical flag indicating wether gene sets tested for discriminating among patients, or for time trends over time.}
\item{ntg}{the number of treatment groups.}
\item{ngs}{the number of tested gene sets.}
\item{nsignif}{the number of significant gene sets at a 5\% FDR (using the default Benjamini & Yekutieli step-up procedure).}
}

\author{
Boris P. Hejblum
}

\seealso{
\code{\link{TcGSA.LR}}
}
\examples{
data(data_simu_TcGSA)

tcgsa_sim_1grp <- TcGSA.LR(expr = expr_1grp, gmt = gmt_sim, Patient_ID = design$Patient_ID,
                           TimePoint = design$TimePoint, func = "linear")
summary(tcgsa_sim_1grp)

tcgsa_sim_2grp <- TcGSA.LR(expr = expr_2grp, gmt = gmt_sim, Patient_ID = design$Patient_ID,
                           TimePoint = design$TimePoint, func = "linear")
summary(tcgsa_sim_2grp)
}

