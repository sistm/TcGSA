% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFit.GS.R
\name{plotFit.GS}
\alias{plotFit.GS}
\title{Plotting function for exploring the fitness of the mixed modeling used in TcGSA}
\usage{
plotFit.GS(
  x,
  expr,
  design,
  subject_name = "Patient_ID",
  time_name = "TimePoint",
  colnames_ID,
  plot_type = c("Fit", "Residuals Obs", "Residuals Est", "Histogram Obs"),
  GeneSetsList,
  color = c("genes", "time", "subjects"),
  marginal_hist = TRUE,
  gg.add = list(theme())
)
}
\arguments{
\item{x}{a \bold{tcgsa} object for \code{clustTrend}, or a
\bold{\link{ClusteredTrends}} object for \code{print.ClusteredTrends} and
\code{plot.ClusteredTrends}.}

\item{expr}{a matrix or dataframe of gene expression.  Its dimension are
\eqn{n}x\eqn{p}, with the \eqn{p} samples in column and the \eqn{n} genes in
row.}

\item{design}{a matrix or dataframe containing the experimental variables that used in the model,
namely \code{subject_name}, \code{time_name}, and \code{covariates_fixed} 
and \code{time_covariates} if applicable.  Its dimension are \eqn{p}x\eqn{m} 
and its row are is in the same order as the columns of \code{expr}.}

\item{subject_name}{the name of the factor variable from \code{design} that contains the information on 
the repetition units used in the mixed model, such as the patient identifiers for instance.  
Default is \code{'Patient_ID'}.  See Details.}

\item{time_name}{the name of a numeric variable from \code{design} that contains 
the information on the time replicates (the time points at which gene 
expression was measured).  Default is \code{'TimePoint'}.  See Details.}

\item{colnames_ID}{the name of the variable from \code{design} that contains the column names of the \code{expr}
expression data matrix.  See Details.}

\item{plot_type}{a character string indicating the type of plot to be drawn.  The options are 
\code{'Fit'}, \code{'Residuals Obs'}, \code{'Residuals Est'} or \code{'Histogram Obs'}.}

\item{GeneSetsList}{a character string containing the names of the gene set whose fit is being checked. 
If several gene sets are being checked, can be a character list or vector of the 
names of those gene sets.}

\item{color}{a character string indicating which color scale should be used. One of the 3 : 
\code{'genes'}, \code{'time'}, \code{'subjects'}, otherwise, no coloring is used.}

\item{marginal_hist}{a logical flag indicating whether marginal histograms should be drawn. 
Only used for \code{'Fit'} plot type. Default is \code{'TRUE'}}

\item{gg.add}{A list of instructions to add to the \code{ggplot2} instructions.
See \code{\link[ggplot2]{+.gg}}.  Default is \code{list(theme())}, which adds nothing
to the plot.}
}
\description{
This function plots graphs informing on the fit of the mixed modeling
of the gene expression performed in TcGSA, for 1 or several gene sets.
}
\examples{

if(interactive()){

data(data_simu_TcGSA)

tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design, 
                          subject_name="Patient_ID", time_name="TimePoint",
                          time_func="linear", crossedRandom=FALSE)
plotFit.GS(x=tcgsa_sim_1grp, expr=expr_1grp, design=design,
				 subject_name="Patient_ID", time_name="TimePoint",
				 colnames_ID="Sample_name", 
				 plot_type="Residuals Obs", 
				 GeneSetsList=c("Gene set 1", "Gene set 2", "Gene set 3",
				                "Gene set 4", "Gene set 5"),
				 color="genes", gg.add=list(guides(color=FALSE))
)

plotFit.GS(x=tcgsa_sim_1grp, expr=expr_1grp, design=design,
          subject_name="Patient_ID", time_name="TimePoint",
          colnames_ID="Sample_name", 
          plot_type="Histogram Obs", 
          GeneSetsList=c("Gene set 1", "Gene set 5"),
          color="genes", gg.add=list(guides(fill=FALSE))
          )
          
plotFit.GS(x=tcgsa_sim_1grp, expr=expr_1grp, design=design,
          subject_name="Patient_ID", time_name="TimePoint",
          colnames_ID="Sample_name", 
          plot_type="Histogram Obs", 
          GeneSetsList=c("Gene set 1", "Gene set 2", "Gene set 3",
		                "Gene set 4", "Gene set 5"),
          color="genes")
}

}
\references{
Hejblum BP, Skinner J, Thiebaut R, (2015) 
Time-Course Gene Set Analysis for Longitudinal Gene Expression Data. 
\emph{PLOS Comput. Biol.} 11(6):e1004310.
doi: 10.1371/journal.pcbi.1004310
}
\seealso{
\code{\link{plot1GS}}, \code{\link{plotSelect.GS}}
}
\author{
Boris P. Hejblum
}
