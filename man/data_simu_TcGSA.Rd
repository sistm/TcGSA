\name{data_simu_TcGSA}
\docType{data}
\alias{data_simu_TcGSA}
\alias{expr_1grp}
\alias{expr_2grp}
\alias{design}
\alias{gmt_sim}

\title{
Simulated Data for TcGSA
}

\description{
Simulated data for 5 gene sets of 50 genes each.  Gene expression is simulated at 5 time points for 10 patients.
}

\usage{
data(data_simu_TcGSA)
}


\value{
  \item{expr_1grp}{
See Details.
}
  \item{expr_2grp}{
See Details.
}
  \item{design}{
a data frame with 5 variables:\itemize{
  \item \code{Patient_ID}: a factor that contains the patient ID.

  \item \code{TimePoint}: a numeric vector or a factor that contains the time points at which gene expression was measured.

\item \code{sample_name}: a character vector with the names of the sample (corresponding to the names of the columns of \code{expr_1grp} and of \code{expr_2grp}).

\item \code{group.var}: a factor that indicates to which of the 2 treatment groups each sample belongs to.

\item \code{Group_paired_ID}
a random paired identifier for paired couples (one in each of the 2 treatment groups) of patients.
}
}
  \item{gmt_sim}{
a \bold{gmt} object containing the gene sets definition.  See \code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}} and \href{https://docs.gsea-msigdb.org/#GSEA/Data_Formats/}{GMT definition on https://docs.gsea-msigdb.org/}.
}
}

\details{

In \code{expr_1grp} all patients belong to the same unique treatment group.  The first 2 gene sets are simulated under the null hypothesis.  The gene sets 3 and 4 are simulated under the alternative hypothesis that there is a significant homogeneous time trend within the gene set.  The gene set 5 is simulated under the alternative hypothesis that there are significant heterogeneous time trends within the gene set.
\cr

In \code{expr_2grp} all patients belong to 2 treatment groups.  The 5 first patients belong to the treatment group '\code{T}', The 5 other patients belong to the treatment group '\code{C}'.  The first 2 gene sets are simulated under the null hypothesis that there is no difference in the time trend between the 2 treatment groups.  The gene sets 3 and 4 are simulated under the alternative hypothesis that there are significantly different homogeneous time trends within the gene set between the 2 treatment groups.  The gene set 5 is simulated under the alternative hypothesis that there are significantly different heterogeneous time trends between the 2 treatment groups within the gene set.

}

\source{
This is simulated data.
}

\author{
Boris P. Hejblum
}

\seealso{
\code{\link{TcGSA.LR}}
}

\examples{

data(data_simu_TcGSA)
summary(expr_1grp)
summary(design)
gmt_sim

}
