#'Plotting Multiple Gene Sets in a single plot
#'
#'This function can plot different representations of the gene expression 
#'in a list of gene sets.
#'
#'If \code{expr} is a matrix or a dataframe, then the "original" data are
#'plotted.  On the other hand, if \code{expr} is a list returned in the
#'\code{'Estimations'} element of \code{\link{TcGSA.LR}}, then it is those
#'"estimations" made by the \code{\link{TcGSA.LR}} function that are plotted.
#'
#'If \code{indiv} is 'genes', then each line of the plot is the median of a
#'gene expression over the patients. On the other hand, if \code{indiv} is
#''patients', then each line of the plot is the median of a patient genes
#'expression in this gene set.
#'
#'This function uses the Gap statistics to determine the optimal number of
#'clusters in the plotted gene set.  See
#'\code{\link[cluster:clusGap]{clusGap}}.
#'
#'@param genesets_list a list of the character strings giving the names of the 
#'genesets to be plotted as they appear in \code{gmt}.
#'
#'@param ncolumns the number of columns used to display the multiple plots. 
#'Default is \code{1}.
#'
#'@param labels List of labels to be added to the plots. 
#'You can also set \code{labels="AUTO"} to auto-generate upper-case labels (such as \code{A}, \code{B}, ...) 
#'or \code{labels="auto"} to auto-generate lower-case labels.
#'Default is \code{NULL}
#'
#'@param expr 
#'either a matrix or dataframe of gene expression upon which
#'dynamics are to be calculated, or a list of gene sets estimation of gene
#'expression.  In the case of a matrix or dataframe, its dimension are \eqn{n}
#'x \eqn{p}, with the \eqn{p} sample in column and the \eqn{n} genes in row.
#'In the case of a list, its length should correspond to the number of gene
#'sets under scrutiny and each element should be an 3 dimension array of
#'estimated gene expression, such as for the list returned in the
#'\code{'Estimations'} element of \code{\link{TcGSA.LR}}.  See details.
#'
#'@param gmt 
#'a \bold{gmt} object containing the gene sets definition.  See
#'\code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}} and
#'definition on \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{www.broadinstitute.org}.
#'
#'@param Subject_ID 
#'a factor of length \eqn{p} that is in the same order as the
#'columns of \code{expr} (when it is a dataframe) and that contains the patient
#'identifier of each sample.
#TODO See Details.
#'
#'@param TimePoint 
#'a numeric vector or a factor of length \eqn{p} that is in
#'the same order as \code{Subject_ID} and the columns of \code{expr} (when it is
#'a dataframe), and that contains the time points at which gene expression was
#'measured.
#TODO See Details.
#'
#'@param baseline 
#'a character string which is the value of \code{TimePoint}
#'that can be used as a baseline.  Default is \code{NULL}, in which case no
#'time point is used as a baseline value for gene expression.  Has to be
#'\code{NULL} when comparing two treatment groups.
#TODO See Details.
#'
#'@param group.var 
#'in the case of several treatment groups, this is a factor of
#'length \eqn{p} that is in the same order as \code{Timepoint},
#'\code{Subject_ID} and the columns of \code{expr}.  It indicates to which
#'treatment group each sample belongs to.  Default is \code{NULL}, which means
#'that there is only one treatment group.  
#TODO See Details.
#'
#'@param Group_ID_paired 
#'a character vector of length \eqn{p} that is in the
#'same order as \code{Timepoint}, \code{Subject_ID}, \code{group.var} and the
#'columns of \code{expr}.  This argument must not be \code{NULL} in the case of
#'a paired analysis, and must be \code{NULL} otherwise.  Default is
#'\code{NULL}.  
#TODO See Details.
#'
#'@param ref 
#'the group which is used as reference in the case of several
#'treatment groups.  Default is \code{NULL}, which means that reference is the
#'first group in alphabetical order of the labels of \code{group.var}.  See
#'Details.
#'
#'@param group_of_interest 
#'the group of interest, for which dynamics are to be
#'computed in the case of several treatment groups.  Default is \code{NULL},
#'which means that group of interest is the second group in alphabetical order
#'of the labels of \code{group.var}.  
#TODO See Details.
#'
#'@param FUNcluster 
#'a function which accepts as first argument a matrix
#'\code{x} and as second argument the number of clusters desired \code{k}, and
#'which returns a list with a component named \code{'cluster'} which is a
#'vector of length \code{n = nrow(x)} of integers in 1:k, determining the clustering
#'or grouping of the n observations.  Default is \code{NULL}, in which case a
#'hierarchical clustering is performed via the function
#'\code{\link[cluster:agnes]{agnes}}, using the metric \code{clustering_metric}
#'and the method \code{clustering_method}.  See \code{'FUNcluster'} in
#'\code{\link[cluster:clusGap]{clusGap}} and Details.
#'
#'@param clustering_metric 
#'character string specifying the metric to be used
#'for calculating dissimilarities between observations in the hierarchical
#'clustering when \code{FUNcluster} is \code{NULL}.  The currently available
#'options are \code{"euclidean"} and \code{"manhattan"}.  Default is
#'\code{"euclidean"}.  See \code{\link[cluster:agnes]{agnes}}.  Also, a \code{"sts"} option 
#'is available in TcGSA.  It implements the 'Short Time Series' distance 
#'[Moller-Levet et al., Fuzzy Clustering of short time series and unevenly distributed 
#'sampling points, \emph{Advances in Intelligent Data Analysis V}:330-340 Springer, 2003]
#'designed specifically for clustering time series.
#'
#'@param clustering_method 
#'character string defining the agglomerative method
#'to be used in the hierarchical clustering when \code{FUNcluster} is
#'\code{NULL}.  The six methods implemented are \code{"average"} ([unweighted
#'pair-]group average method, UPGMA), \code{"single"} (single linkage),
#'\code{"complete"} (complete linkage), \code{"ward"} (Ward's method),
#'\code{"weighted"} (weighted average linkage).  Default is \code{"ward"}.  See
#'\code{\link[cluster:agnes]{agnes}}.
#'
#'@param B 
#'integer specifying the number of Monte Carlo ("bootstrap") samples
#'used to compute the gap statistics.  Default is \code{500}.  See
#'\code{\link[cluster:clusGap]{clusGap}}.
#'
#'@param max_trends 
#'integer specifying the maximum number of different clusters
#'to be tested.  Default is \code{4}.
#'
#'@param aggreg.fun 
#'a character string such as \code{"median"} or \code{"mean"}
#'or the name of any other defined statistics function that returns a single
#'numeric value.  It specifies the function used to aggregate the observations
#'before the clustering. Default is to \code{"median"}.
#'
#'@param trend.fun 
#'a character string such as \code{"mean"} or
#'the name of any other function that returns a single numeric value.  It
#'specifies the function used to calculate the trends of the identified
#'clustered.  Default is to \code{"mean"}.
#'
#'@param methodOptiClust 
#'character string indicating how the "optimal" number
#'of clusters is computed from the gap statistics and their standard
#'deviations. Possible values are \code{"globalmax"}, \code{"firstmax"},
#'\code{"Tibs2001SEmax"}, \code{"firstSEmax"} and \code{"globalSEmax"}.
#'Default is \code{"firstSEmax"}.  See \code{'method'} in
#'\code{\link[cluster:clusGap]{clusGap}}, Details and \emph{Tibshirani et al.,
#'2001} in References.
#'
#'@param indiv a character string indicating by which unit observations are
#'aggregated (through \code{aggreg.fun}) before the clustering.  Possible
#'values are \code{"genes"} or \code{"patients"}.  Default is \code{"genes"}.
#'See Details.
#'
#'@param verbose 
#'logical flag enabling verbose messages to track the computing
#'status of the function.  Default is \code{TRUE}.
#'
#'@param clustering 
#'logical flag.  If \code{FALSE}, there is no clustering
#'representation; if \code{TRUE}, the lines are colored according to which
#'cluster they belong to.  Default is \code{TRUE}.  See Details.
#'
#'@param showTrend 
#'logical flag.  If \code{TRUE}, a black line is added for
#'each cluster, representing the corresponding \code{trend.fun}.  Default is
#'\code{TRUE}.
#'
#'@param smooth 
#'logical flag.  If \code{TRUE} and \code{showTrend} is also
#'\code{TRUE}, the representation of each cluster \code{trend.fun} is smoothed
#'using cubic polynomials (see \code{\link[ggplot2:geom_smooth]{geom_smooth}}.
#'Default is \code{TRUE}. 
#'At the moment, must accept parameter \code{"na.rm"} (which is automatically set to \code{TRUE}). 
#'This might change in future versions
#'
#'@param time_unit 
#'the time unit to be displayed (such as \code{"Y"},
#'\code{"M"}, \code{"W"}, \code{"D"}, \code{"H"}, etc) next to the values of
#'\code{TimePoint} on the x-axis.  Default is \code{""}, in which case the time 
#'scale on the x-axis is proportional to the time values.
#'
#'@param y.lab 
#'character specifying the annotation of the y axis.  If \code{NULL}, an
#'annotation is automatically generated, if \code{""}, no annotation appears.  Default is
#'\code{NULL}.
#'
#'@param desc 
#'a logical flag. If \code{TRUE}, a line is added to the title of
#'the plot with the description of the gene set plotted (from the gmt file).
#'Default is \code{TRUE}.
#'
#'@param lab.cex 
#'a numerical value giving the amount by which lab labels text
#'should be magnified relative to the default \code{1}.
#'
#'@param axis.cex 
#'a numerical value giving the amount by which axis annotation
#'text should be magnified relative to the default \code{1}.
#'
#'@param main.cex 
#'a numerical value giving the amount by which title text
#'should be magnified relative to the default \code{1}.
#'
#'@param margins
#'a numerical value giving the amount by which the margins
#'should be reduced or increased relative to the default \code{1}.
#'
#'@param line.size
#'a numerical value giving the amount by which the line sizes
#'should be reduced or increased relative to the default \code{1}.
#'
#'@param y.lab.angle 
#'a numerical value (in [0, 360]) giving the orientation by
#'which y-label text should be turned (anti-clockwise).  Default is \code{90}.
#'See \code{\link{element_text}}.
#'
#'@param x.axis.angle 
#'a numerical value (in [0, 360]) giving the orientation by
#'which x-axis annotation text should be turned (anti-clockwise).  Default is
#'\code{45}.
#'
#'@param y.lim 
#'a numeric vector of length 2 giving the range of the y-axis.
#'See \code{\link{plot.default}}.
#'
#'@param x.lim 
#'if numeric, will create a continuous scale, if factor or
#'character, will create a discrete scale.  Observations not in this range will
#'be dropped.  See \code{\link{xlim}}.
#'
#'@param gg.add 
#'A list of instructions to add to the \code{ggplot2} instructions.  
#'See \code{\link{+.gg}}.  Default is \code{list(theme())}, which adds nothing
#'to the plot.
#'
#'@param show_plot 
#'logical flag. If \code{FALSE}, no plot is drawn. Default is \code{TRUE}.
#'
#'@return A list with 2 elements:\itemize{
#'   \item \code{classif}: a \code{data.frame} with  the 2 following variables: \code{ProbeID} which 
#'   contains the IDs of the probes of the plotted gene set, and \code{Cluster} containing $
#'   which cluster the probe belongs to. If \code{clustering} is \code{FALSE}, then \code{Cluster} is \code{NA} for all the probes.
#'   \item \code{p}: a \code{ggplot} object containing the plot
#'}
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link[ggplot2:ggplot]{ggplot}}, \code{\link[cluster:clusGap]{clusGap}}
#'
#'@import ggplot2
#'
#'@importFrom cowplot plot_grid
#'
#'@export
#'
#'
plotMultipleGS <- function(genesets_list, ncolumns=1, labels=NULL,
						   expr, gmt, Subject_ID, TimePoint, 
						   baseline=NULL,
						   group.var=NULL, Group_ID_paired=NULL, ref=NULL, group_of_interest=NULL,
						   FUNcluster=NULL, clustering_metric="euclidian", clustering_method="ward", B=500,
						   max_trends=4, aggreg.fun="median", trend.fun="median",
						   methodOptiClust = "firstSEmax",
						   indiv="genes",
						   verbose=TRUE,
						   clustering=TRUE, showTrend=TRUE, smooth=TRUE,
						   time_unit="", y.lab=NULL, desc=TRUE,
						   lab.cex=1, axis.cex=1, main.cex=1, y.lab.angle=90, x.axis.angle=45, margins=1, line.size=1,
						   y.lim=NULL, x.lim=NULL, 
						   gg.add=list(),
						   show_plot=TRUE){
	
	plot_list <- list()
	for(gs in genesets_list){
		id_gs <- which(gmt$geneset.names == gs)
		if(desc){
			mytitle <- paste0(gmt$geneset.names[id_gs], ": ", gmt$geneset.descriptions[id_gs])#, " (", gmt$geneset.descriptions[id_gs], ")")
		}else{
			mytitle <- gmt$geneset.names[id_gs]
		}
		plot_list[[gs]] <- plot1GS(expr=expr, gmt=gmt,
								   Subject_ID = Subject_ID, TimePoint = TimePoint,
								   geneset.name = gs,
								   baseline = baseline, group.var = group.var, Group_ID_paired = Group_ID_paired,
								   ref = ref, group_of_interest = group_of_interest,
								   FUNcluster = FUNcluster, clustering_metric = clustering_metric, 
								   clustering_method=clustering_method, B=B,
								   max_trends = max_trends, 
								   aggreg.fun = aggreg.fun, trend.fun = trend.fun,
								   methodOptiClust = methodOptiClust,
								   indiv = indiv, verbose=verbose, 
								   clustering = clustering, showTrend = showTrend, smooth = smooth,
								   time_unit = time_unit,
								   y.lab=y.lab, desc=desc,
								   lab.cex = lab.cex, axis.cex = axis.cex, main.cex = main.cex, 
								   y.lab.angle=y.lab.angle, x.axis.angle=x.axis.angle, margins=margins, line.size=line.size,
								   y.lim=y.lim, x.lim=x.lim, 
								   title = mytitle,
								   gg.add = gg.add,
								   plot=FALSE
		)$p
	}
	if(show_plot){
		print(cowplot::plot_grid(plotlist=plot_list, ncol=ncolumns), labels = labels)
	}
	invisible(plot_list)
}


