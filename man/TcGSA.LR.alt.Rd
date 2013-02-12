\name{TcGSA.LR.alt}
\alias{TcGSA.LR.alt}

\title{
Computing the Likelihood Ratios for the Gene Sets under Scrutiny
}
\description{
This function computes the Likelihood Ratios for the gene sets under scrutiny, as well as estimations of genes dynamics inside those gene sets through mixed models.
}
\usage{
TcGSA.LR.alt(expr, gmt, Patient_ID, TimePoint, 
         func = "linear",  minGSsize=10, maxGSsize=500, 
         group.var = NULL, separatePatients=FALSE)
}

\arguments{
  \item{expr}{
a matrix or dataframe of gene expression.  Its dimension are \eqn{n}x\eqn{p}, with the \eqn{p} samples in column and the \eqn{n} genes in row.
}
  \item{gmt}{
a \bold{gmt} object containing the gene sets definition.  See \code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}} and \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}{GMT definition on www.broadinstitute.org}.
}
  \item{Patient_ID}{
a factor of length \eqn{p} that is in the same order as the columns of \code{expr} (when it is a dataframe) and that contains the patient identifier of each sample.%  See Details.
}
  \item{TimePoint}{
a numeric vector or a factor of length \eqn{p} that is in the same order as \code{TimePoint} and the columns of \code{expr} (when it is a dataframe), and that contains the time points at which gene expression was measured.%  See Details.
}
  \item{func}{
the form of the time trend. Can be either one of \code{"linear"}, \code{"cubic"} or \code{"splines"}. The \code{"splines"} form corresponds to the natural cubic B-splines (see also \code{\link[splines:ns]{ns}}).  If there are only a few timepoints, a \code{"linear"} form should be sufficient. Otherwise, the \code{"cubic"} form is more parsimonious than the \code{"splines"} form, and should be sufficiently flexible.
}  
\item{minGSsize}{
the minimum number of genes in a gene set.  If there are less genes than this number in one of the gene sets under scrutinity, the Likelihood Ratio of this gene set is not computed (the mixed model are not fitted). Default is \code{10} genes as the minimum.
}
  \item{maxGSsize}{
the maximum number of genes in a gene set.  If there are more genes than this number in one of the gene sets under scrutinity, the Likelihood Ratio of this gene set is not computed (the mixed model are not fitted).  This is to avoid very long computation times.  Default is \code{500} genes as the maximum.
}
  \item{group.var}{
in the case of several treatment groups, this is a factor of length \eqn{p} that is in the same order as \code{TimePoint}, \code{Patient_ID}, \code{sample_name} and the columns of \code{expr}.  It indicates to which treatment group each sample belongs to.  Default is \code{NULL}, which means that there is only one treatment group.  See Details.
}
  \item{separatePatients}{
logical flag indicating that the analysis identifies gene sets that discriminates patients rather than gene sets than have a significant trend over time.  Default is \code{FALSE}.  See Details.
}
}

\details{
This Time-course Gene Set Analysis aims at identifying gene sets that are not stable over time, either homogeneously  or heterogeneously (see \emph{Hejblum et al, 2012})in terms of their probes.  And when the argument \code{separatePatients} is \code{TRUE}, instead of identifying gene sets that have a significant trend over time, \emph{TcGSA} identifies gene sets that have significantly different trends over time depending on the patient.

%TODO
}

\value{
\code{TcGSA.LR} returns a \code{tcgsa} object, which is a list with the 5 following elements:
  \item{fit}{a data frame that contains the 3 following variables:
    \itemize{
      \item \code{LR}: the likelihood ratio between the model under the null hypothesis and the model under the alternative hypothesis.
      \item \code{CVG_H0}: convergence status of the model under the null hypothesis.
      \item \code{CVG_H1}: convergence status of the model under the alternative hypothesis.
    }
  }
  \item{func_form}{
    a character string passing along the value of the \code{func} argument used in the call.
  }
  \item{GeneSets_gmt}{
    a \code{gmt} object passing along the value of the \code{gmt} argument used in the call.
  }
  \item{group.var}{
    a factor passing along the value of the \code{groupe.var} argument used in the call.
  }
  \item{separatePatients}{
    a logical flag passing along the value of the \code{separatePatients} argument used in the call.
  }
  \item{Estimations}{
    a list of 3 dimensions arrays.  Each element of the list (i.e. each array) corresponds to the estimations of gene expression dynamics for each of the gene sets under scrutiny (obtained from mixed models).  The first dimension of those arrays is the genes included in the concerned gene set, the second dimension is the \code{Patient_ID}, and the third dimension is the \code{TimePoint}.  The values inside those arrays are estimated gene expressions..
  }
  \item{splines_DF}{
  	
  }
}

\references{
Hejblum, B.P., Skinner, J., Thiebaut, R., 2012, TcGSA: Time-course gene set analysis of longitudinal gene expression, with an application to a therapeutic HIV vaccine trial, \bold{submitted}.
}

\author{
Boris P. Hejblum
}

\seealso{
  \code{\link{summary.TcGSA}}, \code{\link{plot.TcGSA}}
}
\examples{
data(data_simu_TcGSA)

tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, Patient_ID=design$Patient_ID,
                           TimePoint=design$TimePoint, func="linear")
tcgsa_sim_1grp
summary(tcgsa_sim_1grp)
plot(x=tcgsa_sim_1grp, expr=expr_1grp, 
     Patient_ID=design$Patient_ID, TimePoint=design$TimePoint,
     baseline=1, 
     B=100,
     time_unit="H"
     )
     
     
tcgsa_sim_2grp <- TcGSA.LR(expr=expr_2grp, gmt=gmt_sim, Patient_ID=design$Patient_ID,
                           TimePoint=design$TimePoint, func="linear", 
                           group.var=design$group.var)
tcgsa_sim_2grp
}
