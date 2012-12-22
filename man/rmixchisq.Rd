\name{rmixchisq}
\alias{rmixchisq}

\title{
Random Generation of Chi-square Mixtures
}

\description{
\code{rmixchisq} is used to simulate a mixture of chi-square distributions that corresponds to the null distribution of the Likelihood Ratio between 2 nested mixed models.
}

\usage{
rmixchisq(n, s, q)
}

\arguments{
  \item{n}{
        number of observations.
}
  \item{s}{
        number of fixed effects to be tested.
}
  \item{q}{
        number of random effects to be tested.
}
}

\details{
The approximate null distribution of a likelihood ratio for 2 nested mixed models, where both fixed and random effects are tested simultaneously, is a very specific mixture of \eqn{\chi^2} distributions [\cite{Self & Liang (1987), Stram & Lee (1994) and Stram & Lee (1995)}].  It depends on both the number of random effects and the number of fixed effects to be tested simultaneously.
}
\value{
A vector of random independent observations of the chisquare mixture identified by the values of \code{s} and \code{q}.
}
\references{
Self, S. G. and Liang, K., 1987, Asymptotic properties of maximum likelihood estimators and likelihood ratio tests under nonstandard conditions, \emph{Journal of the American Statistical Association} 82: 605--610.

Stram, D. O. and Lee, J. W., 1994, Variance components testing in the longitudinal mixed effects model, \emph{Biometrics} 50: 1171--1177.

Stram, D. O. and Lee, J. W., 1995, Corrections to "Variance components testing in the longitudinal mixed effects model" by Stram, D. O. and Lee, J. W.; 50: 1171--1177 (1994), \emph{Biometrics} 51: 1196.
}

\author{
Boris P. Hejblum
}

\seealso{
\code{\link{pval_simu}}
}

\examples{
sample_mixt <- rmixchisq(n=1000, s=3, q=3)
plot(density(sample_mixt))
}
