\name{pval_simu}
\alias{pval_simu}
\title{Computing P-values with a Simulated Sample from the Null Distribution}
\usage{
  pval_simu(s, theo_dist)
}
\arguments{
  \item{s}{the observation whose p-value is computed.  For
  instance a Likelihood Ratio.}

  \item{theo_dist}{the sample of the distribution under the
  null hypothesis.}
}
\value{
  The p-value associated to the observation \code{s}.
}
\description{
  This function computes the p-value of a statistic using a
  simulated sample from its therotical null distribution.
}
\examples{
theo_dist <- rnorm(n=10000, mean=0, sd=1)
pval_simu(s=1.96, theo_dist)
1-pnorm(q=1.96, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
}
\author{
  Boris P. Hejblum
}
\seealso{
  \code{\link{rmixchisq}}
}

