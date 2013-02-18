\name{TcGSA.LR}
\alias{print.TcGSA}
\alias{print.TcGSA.LR}
\alias{TcGSA.LR}
\title{Computing the Likelihood Ratios for the Gene Sets under Scrutiny}
\usage{
  TcGSA.LR(expr, gmt, design, subject_name = "Patient_ID",
    time_name = "TimePoint", crossedRandom = FALSE,
    covariates_fixed = "", time_covariates = "",
    time_func = "linear", group_name = "",
    separateSubjects = FALSE, minGSsize = 10,
    maxGSsize = 500)

  \method{print}{TcGSA} (x, ...)
}
\arguments{
  \item{expr}{a matrix or dataframe of gene expression.
  Its dimension are \eqn{n}x\eqn{p}, with the \eqn{p}
  samples in column and the \eqn{n} genes in row.}

  \item{gmt}{a \bold{gmt} object containing the gene sets
  definition.  See
  \code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}} and
  definition on
  \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29GMT}{www.broadinstitute.org}.}

  \item{design}{a matrix or dataframe containing the
  experimental variables that used in the model, namely
  \code{subject_name}, \code{time_name}, and
  \code{covariates_fixed} and \code{time_covariates} if
  applicable.  Its dimension are \eqn{p}x\eqn{m} and its
  row are is in the same order as the columns of
  \code{expr}.}

  \item{subject_name}{the name of the factor variable from
  \code{design} that contains the information on the
  repetition units used in the mixed model, such as the
  patient identifiers for instance. Default is
  \code{'Patient_ID'}.  See Details.}

  \item{time_name}{the name of the numeric or factor
  variable from \code{design} contains the information on
  the time replicates (the time points at which gene
  expression was measured).  Default is \code{'TimePoint'}.
  See Details.}

  \item{crossedRandom}{logical flag indicating wether the
  random effects of the subjects and of the time points
  should be modeled as one crossed random effect or as two
  separated random effects. Default is \code{FALSE}. See
  details.}

  \item{covariates_fixed}{a character vector with the names
  of numeric or factor variables from the \code{design}
  matrix that should appear as fixed effects in the model.
  See details. Default is \code{""}, which corresponds to
  no covariates in the model.}

  \item{time_covariates}{a character vector with the names
  of numeric or factor variables from the \code{design}
  matrix that should appear as fixed effects interaction
  with the \code{time_name} variable in the model.  See
  details.  Default is \code{""}, which corresponds to no
  covariates in the model.}

  \item{time_func}{the form of the time trend. Can be
  either one of \code{"linear"}, \code{"cubic"},
  \code{"splines"} or specified by the user as an
  expression using names of variables from the
  \code{design} matrix. The \code{"splines"} form
  corresponds to the natural cubic B-splines (see also
  \code{\link[splines:ns]{ns}}).  If there are only a few
  timepoints, a \code{"linear"} form should be sufficient.
  Otherwise, the \code{"cubic"} form is more parsimonious
  than the \code{"splines"} form, and should be
  sufficiently flexible.}

  \item{minGSsize}{the minimum number of genes in a gene
  set.  If there are less genes than this number in one of
  the gene sets under scrutinity, the Likelihood Ratio of
  this gene set is not computed (the mixed model are not
  fitted). Default is \code{10} genes as the minimum.}

  \item{maxGSsize}{the maximum number of genes in a gene
  set.  If there are more genes than this number in one of
  the gene sets under scrutinity, the Likelihood Ratio of
  this gene set is not computed (the mixed model are not
  fitted).  This is to avoid very long computation times.
  Default is \code{500} genes as the maximum.}

  \item{group_name}{in the case of several treatment
  groups, the name of a factor variable from the
  \code{design} matrix.  It indicates to which treatment
  group each sample belongs to.  Default is \code{""},
  which means that there is only one treatment group.  See
  Details.}

  \item{separateSubjects}{logical flag indicating that the
  analysis identifies gene sets that discriminates patients
  rather than gene sets than have a significant trend over
  time.  Default is \code{FALSE}.  See Details.}

  \item{x}{an object of class '\code{TcGSA}'.}

  \item{...}{further arguments passed to or from other
  methods.}
}
\value{
  \code{TcGSA.LR} returns a \code{tcgsa} object, which is a
  list with the 5 following elements: \itemize{ \item fit a
  data frame that contains the 3 following variables:
  \itemize{ \item \code{LR}: the likelihood ratio between
  the model under the null hypothesis and the model under
  the alternative hypothesis. \item \code{CVG_H0}:
  convergence status of the model under the null
  hypothesis. \item \code{CVG_H1}: convergence status of
  the model under the alternative hypothesis. } \item
  \code{time_func}: a character string passing along the
  value of the \code{func} argument used in the call. \item
  \code{GeneSets_gmt}: a \code{gmt} object passing along
  the value of the \code{gmt} argument used in the call.
  \item \code{group.var}: a factor passing along the
  \code{group_name} variable from the \code{design} matrix.
  \item \code{separateSubjects}: a logical flag passing
  along the value of the \code{separateSubjects} argument
  used in the call. \item \code{Estimations}: a list of 3
  dimensions arrays.  Each element of the list (i.e. each
  array) corresponds to the estimations of gene expression
  dynamics for each of the gene sets under scrutiny
  (obtained from mixed models).  The first dimension of
  those arrays is the genes included in the concerned gene
  set, the second dimension is the \code{Patient_ID}, and
  the third dimension is the \code{TimePoint}.  The values
  inside those arrays are estimated gene expressions. \item
  \code{time_DF}: the degree of freedom of the natural
  splines functions }
}
\description{
  This function computes the Likelihood Ratios for the gene
  sets under scrutiny, as well as estimations of genes
  dynamics inside those gene sets through mixed models.
}
\details{
  This Time-course Gene Set Analysis aims at identifying
  gene sets that are not stable over time, either
  homogeneously or heterogeneously (see \emph{Hejblum et
  al, 2012})in terms of their probes.  And when the
  argument \code{separateSubjects} is \code{TRUE}, instead
  of identifying gene sets that have a significant trend
  over time, \emph{TcGSA} identifies gene sets that have
  significantly different trends over time depending on the
  patient.
}
\examples{
data(data_simu_TcGSA)

tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design,
													 subject_name="Patient_ID", time_name="TimePoint",
                          time_func="linear", crossedRandom=FALSE)
tcgsa_sim_1grp
summary(tcgsa_sim_1grp)
plot(x=tcgsa_sim_1grp, expr=expr_1grp,
    Subject_ID=design$Patient_ID, TimePoint=design$TimePoint,
    baseline=1,
    B=100,
    time_unit="H"
    )


tcgsa_sim_2grp <- TcGSA.LR(expr=expr_2grp, gmt=gmt_sim, design=design,
													 subject_name="Patient_ID", time_name="TimePoint",
                          time_func="linear", crossedRandom=FALSE,
                          group_name="group.var")
tcgsa_sim_2grp
}
\author{
  Boris P. Hejblum
}
\references{
  Hejblum, B.P., Skinner, J., Thiebaut, R., 2012, TcGSA:
  Time-course gene set analysis of longitudinal gene
  expression, with an application to a therapeutic HIV
  vaccine trial, \bold{submitted}.
}
\seealso{
  \code{\link{summary.TcGSA}}, \code{\link{plot.TcGSA}}
}

