#'Extracting a subset of a TcGSA object
#'
#'This function extracts a subset of a TcGSA object from a list of gene sets
#'
#'@param tcgsa 
#'a \code{tcgsa} object.
#'
#'@param list_tcgsa_subset 
#'a vector of the gene set names for the subset
#'
#'@return \code{TcGSA.subset} returns a \code{tcgsa} object, which is a list with
#'the 5 following elements:
#'\itemize{
#'\item fit a data frame that contains the 7 following variables:
#'\itemize{ 
#'\item \code{LR}: the likelihood ratio between the model under the
#'null hypothesis and the model under the alternative hypothesis.  
#'\item \code{AIC_H0}: AIC criterion for the model under the null hypothesis.
#'\item \code{AIC_H1}: AIC criterion for the model under the alternative hypothesis.
#'\item \code{BIC_H0}: BIC criterion for the model under the null hypothesis.
#'\item \code{BIC_H1}: BIC criterion for the model the alternative hypothesis.
#'\item \code{CVG_H0}: convergence status of the model under the null hypothesis.
#'\item \code{CVG_H1}: convergence status of the model under the alternative
#'hypothesis.
#'}
#'\item \code{time_func}: a character string passing along the value of the
#'\code{time_func} argument used in the call.
#'\item \code{GeneSets_gmt}: a \code{gmt} object passing along the value of the
#'\code{gmt} argument used in the call.
#'\item \code{group.var}: a factor passing along the \code{group_name} variable
#'from the \code{design} matrix.
#'\item \code{separateSubjects}: a logical flag passing along the value of the
#'\code{separateSubjects} argument used in the call.
#'\item \code{Estimations}: a list of 3 dimensions arrays.  Each element of the
#'list (i.e. each array) corresponds to the estimations of gene expression
#'dynamics for each of the gene sets under scrutiny (obtained from mixed
#'models).  The first dimension of those arrays is the genes included in the
#'concerned gene set, the second dimension is the \code{Patient_ID}, and the
#'third dimension is the \code{TimePoint}.  The values inside those arrays are
#'estimated gene expressions.
#'\item \code{time_DF}: the degree of freedom of the natural splines functions
#'}
#'
#'@author Boris P. Hejblum, Damien Chimits
#'
#'@seealso \code{\link{summary.TcGSA}}, \code{\link{plot.TcGSA}}, 
#'and \code{\link{TcGSA.LR.parallel}} for an implementation using 
#'parallel computing
#'
#'@references Hejblum, B.P., Skinner, J., Thiebaut, R., 2014, TcGSA: a gene set approach for longitudinal gene expression data analysis, \bold{submitted}.
#'
#'@examples
#'
#'data(data_simu_TcGSA)
#'
#'tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design, 
#'                           subject_name="Patient_ID", time_name="TimePoint",
#'                           time_func="linear", crossedRandom=FALSE)
#'                           
#'tcgsa_1gr_subset <- TcGSA.subset(tcgsa_sim_1grp,c("Gene set 3", "Gene set 5"))
#'

TcGSA.subset <- 
	function(tcgsa, list_tcgsa_subset){
		index_subset <- which(tcgsa$GeneSets_gmt$geneset.names %in% list_tcgsa_subset)
		tcgsa_subset <- tcgsa
		tcgsa_subset$fit <- tcgsa_subset$fit[index_subset,]
		tcgsa_subset$GeneSets_gmt$genesets <- tcgsa_subset$GeneSets_gmt$genesets[index_subset]
		tcgsa_subset$GeneSets_gmt$geneset.names <- tcgsa_subset$GeneSets_gmt$geneset.names[index_subset]
		tcgsa_subset$GeneSets_gmt$geneset.descriptions <- tcgsa_subset$GeneSets_gmt$geneset.descriptions[index_subset]
		tcgsa_subset$Estimations <- tcgsa_subset$Estimations[index_subset]
		return(tcgsa_subset)
	}