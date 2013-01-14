\name{TcGSA.LR.parallel}
\alias{TcGSA.LR.parallel}

\title{
Parallel computing the Likelihood Ratios for the Gene Sets under Scrutiny
}

\description{
A parallel versioon of the function \code{\link{TcGSA.LR}} to be used on a cluster of computing processors.  This function computes the Likelihood Ratios for the gene sets under scrutiny, as well as estimations of genes dynamics inside those gene sets through mixed models.
}

\usage{
TcGSA.LR.parallel(Nproc, type_connec, 
                  expr, gmt, Patient_ID, TimePoint,
                  func = "linear", maxGSsize=500, 
                  group.var = NULL, separatePatients=FALSE,
                  monitorfile="")
}

\arguments{
  \item{Nproc}{
The number of processors available on the cluster.
  }
  \item{type_connec}{
The type of connection between the processors. Supported cluster types are \code{"SOCK"}, \code{"PVM"}, \code{"MPI"}, and \code{"NWS"}. See also \code{\link{makeCluster}}.
}
  \item{expr}{
a matrix or dataframe of gene expression. Its dimensions are \eqn{n}x\eqn{p}, with the \eqn{p} sample in column and the \eqn{n} genes in row.
  }
  \item{gmt}{
a \code{gmt} object containing the gene sets definition. See \code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}} and \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}{GMT definition on www.broadinstitute.org}.
  }
  \item{Patient_ID}{
a factor of length \eqn{p} hat is in the same order as the columns of \code{expr} and that contains the patient identifier of each sample.% See Details.
  }
  \item{TimePoint}{
a numeric vector or a factor of length \eqn{p} that is in the same order as \code{TimePoint} and the columns of \code{expr}, and that contains the time points at which gene expression was measured.% See Details.
  }
  \item{func}{
the form of the time trend. Can be either one of \code{"linear"}, \code{"cubic"} or \code{"splines"}. The \code{"splines"} form corresponds to the natural cubic B-splines (see also \code{\link[splines:ns]{ns}}).  If there are only a few timepoints, a \code{"linear"} form should be sufficient.  Otherwise, the \code{"cubic"} form is more parsimonious than the \code{"splines"} form, and should be sufficiently flexible.
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
    \item{monitorfile}{
a writable \link{connection} or a character string naming a file to write into, to monitor the progress of the analysis.  Default is \code{""} which is no monitoring.  See Details.
  }
}

\details{
This Time-course Gene Set Analysis aims at identifying gene sets that are not stable over time, either homogeneously  or heterogeneously (see \emph{Hejblum et al, 2012}) in terms of their probes.  And when the argument \code{separatePatients} is \code{TRUE}, instead of identifying gene sets that have a significant trend over time (possibly with probes heterogeneity of this trend), \emph{TcGSA} identifies gene sets that have significantly different trends over time depending on the patient.

If the \code{monitorfile} argument is a character string naming a file to write into, in the case of a new file that does not exist yet, such a new file will be created. A line is written each time one of the gene sets under scrutiny has been analysed (i.e. the two mixed models have been fitted, see \code{\link{TcGSA.LR}}) by one of the parallelized processors.

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
    a list that contains the 2 following elements:
    \itemize{
      \item \code{FixEf}: a list of the fixed effects estimations, with one element for each gene set modelled.
      \item \code{RanEf}: a list of the random effects estimations, with one element for each gene set modelled.
    }See also \code{\link{estim.TcGSA}}.
  }
}

\references{
Hejblum, B.P., Skinner, J., Thiebaut, R., 2012, TcGSA: Time-course gene set analysis of longitudinal gene expression, with an application to a therapeutic HIV vaccine trial, \bold{submitted}.
}

\author{
Boris P. Hejblum
}

\seealso{
  \code{\link{TcGSA.LR}}, \code{\link{summary.TcGSA}}, \code{\link{estim.TcGSA}}, \code{\link{makeCluster}}
}
\examples{
data(data_simu_TcGSA)

\dontrun{
tcgsa_sim_1grp <- TcGSA.LR.parallel(Nproc = 2, type_connec = 'SOCK', 
                                    expr = expr_1grp, gmt = gmt_sim, 
                                    Patient_ID = design$Patient_ID,  
                                    TimePoint = design$TimePoint, 
                                    func = "linear") # monitorfile="monitorEx.txt"
tcgsa_sim_1grp
}

}
