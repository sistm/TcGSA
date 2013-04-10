\name{plotPat.TcGSA}
\alias{plotPat.TcGSA}
\title{Plot a Gene Set Trends Heatmap for each Patient.}
\usage{
  plotPat.TcGSA(x, threshold = 0.05, myproc = "BY",
    nbsimu_pval = 1e+06, expr, Subject_ID, TimePoint,
    baseline = NULL, only.signif = TRUE, group.var = NULL,
    Group_ID_paired = NULL, ref = NULL,
    group_of_interest = NULL, FUNcluster = NULL,
    clustering_metric = "euclidian",
    clustering_method = "ward", B = 500, max_trends = 4,
    aggreg.fun = "median", methodOptiClust = "firstSEmax",
    verbose = TRUE, clust_trends = NULL, N_clusters = NULL,
    myclusters = NULL, label.clusters = NULL,
    prev_rowCL = NULL, descript = TRUE, plotAll = TRUE,
    color.vec = c("darkred", "#D73027", "#FC8D59", "snow", "#91BFDB", "#4575B4", "darkblue"),
    legend.breaks = NULL, label.column = NULL,
    time_unit = "", cex.label.row = 1,
    cex.label.column = 1, margins = c(5, 25),
    heatKey.size = 1, dendrogram.size = 1,
    heatmap.height = 1, heatmap.width = 1,
    cex.clusterKey = 1, cex.main = 1,
    horiz.clusterKey = TRUE, main = NULL, subtitle = NULL,
    ...)
}
\arguments{
  \item{x}{a \bold{tcgsa} object.}

  \item{threshold}{the threshold at which the FDR or the
  FWER should be controlled.}

  \item{myproc}{a vector of character strings containing
  the names of the multiple testing procedures for which
  adjusted p-values are to be computed. This vector should
  include any of the following: "\code{Bonferroni}",
  "\code{Holm}", "\code{Hochberg}", "\code{SidakSS}",
  "\code{SidakSD}", "\code{BH}", "\code{BY}", "\code{ABH}",
  "\code{TSBH}". See
  \code{\link[multtest:mt.rawp2adjp]{mt.rawp2adjp}} for
  details.  Default is "\code{BY}", the Benjamini &
  Yekutieli (2001) step-up FDR-controlling procedure
  (general dependency structures).  In order to control the
  FWER(in case of an analysis that is more a hypothesis
  confirmation than an exploration of the expression data),
  we recommand to use "\code{Holm}", the Holm (1979)
  step-down adjusted p-values for strong control of the
  FWER.}

  \item{nbsimu_pval}{the number of observations under the
  null distribution to be generated in order to compute the
  p-values. Default is \code{1e+06}.}

  \item{expr}{either a matrix or dataframe of gene
  expression upon which dynamics are to be calculated, or a
  list of gene sets estimation of gene expression.  In the
  case of a matrix or dataframe, its dimension are \eqn{n}
  x \eqn{p}, with the \eqn{p} sample in column and the
  \eqn{n} genes in row. In the case of a list, its length
  should correspond to the number of gene sets under
  scrutiny and each element should be an 3 dimension array
  of estimated gene expression, such as for the list
  returned in the \code{'Estimations'} element of
  \code{\link{TcGSA.LR}}.  See details.}

  \item{Subject_ID}{a factor of length \eqn{p} that is in
  the same order as the columns of \code{expr} (when it is
  a dataframe) and that contains the patient identifier of
  each sample.}

  \item{TimePoint}{a numeric vector or a factor of length
  \eqn{p} that is in the same order as \code{Subject_ID}
  and the columns of \code{expr} (when it is a dataframe),
  and that contains the time points at which gene
  expression was measured.}

  \item{baseline}{a character string which is the value of
  \code{TimePoint} used as baseline.  See Details.}

  \item{only.signif}{logical flag for plotting only the
  significant gene sets. If \code{FALSE}, all the gene sets
  from the \bold{gmt} object contained in \code{x} are
  plotted.  Default is \code{TRUE}.}

  \item{group.var}{in the case of several treatment groups,
  this is a factor of length \eqn{p} that is in the same
  order as \code{Timepoint}, \code{Subject_ID},
  \code{sample_name} and the columns of \code{expr}.  It
  indicates to which treatment group each sample belongs
  to.  Default is \code{NULL}, which means that there is
  only one treatment group.  See Details.}

  \item{Group_ID_paired}{a character vector of length
  \eqn{p} that is in the same order as \code{Timepoint},
  \code{Subject_ID}, \code{sample_name}, \code{group.var}
  and the columns of \code{expr}.  This argument must not
  be \code{NULL} in the case of a paired analysis, and must
  be \code{NULL} otherwise.  Default is \code{NULL}.  See
  Details.}

  \item{ref}{the group which is used as reference in the
  case of several treatment groups.  Default is
  \code{NULL}, which means that reference is the first
  group in alphabetical order of the labels of
  \code{group.var}.  See Details.}

  \item{group_of_interest}{the group of interest, for which
  dynamics are to be computed in the case of several
  treatment groups.  Default is \code{NULL}, which means
  that group of interest is the second group in
  alphabetical order of the labels of \code{group.var}.
  See Details.}

  \item{FUNcluster}{the clustering function used to
  agglomerate genes in trends.  Default is \code{NULL}, in
  which a hierachical clustering is performed via the
  function \code{\link[cluster:agnes]{agnes}}, using the
  metric \code{clustering_metric} and the method
  \code{clustering_method}.  See
  \code{\link[cluster:clusGap]{clusGap}}}

  \item{clustering_metric}{character string specifying the
  metric to be used for calculating dissimilarities between
  observations in the hierarchical clustering when
  \code{FUNcluster} is \code{NULL}.  The currently
  available options are \code{"euclidean"} and
  \code{"manhattan"}.  Default is \code{"euclidean"}.  See
  \code{\link[cluster:agnes]{agnes}}.}

  \item{clustering_method}{character string defining the
  agglomerative method to be used in the hierarchical
  clustering when \code{FUNcluster} is \code{NULL}.  The
  six methods implemented are \code{"average"} ([unweighted
  pair-]group average method, UPGMA), \code{"single"}
  (single linkage), \code{"complete"} (complete linkage),
  \code{"ward"} (Ward's method), \code{"weighted"}
  (weighted average linkage).  Default is \code{"ward"}.
  See \code{\link[cluster:agnes]{agnes}}.}

  \item{B}{integer specifying the number of Monte Carlo
  ("bootstrap") samples used to compute the gap statistics.
  Default is \code{500}.  See
  \code{\link[cluster:clusGap]{clusGap}}.}

  \item{max_trends}{integer specifying the maximum number
  of different clusters to be tested.  Default is
  \code{4}.}

  \item{aggreg.fun}{a character string such as
  \code{"mean"}, \code{"median"} or the name of any other
  statistics function defined that returns a single numeric
  value.  It specifies the function used to aggregate the
  observations before the clustering.  Default is to
  \code{median}.  Default is \code{"median"}.}

  \item{methodOptiClust}{character string indicating how
  the "optimal"" number of clusters is computed from the
  gap statistics and their standard deviations. Possible
  values are \code{"globalmax"}, \code{"firstmax"},
  \code{"Tibs2001SEmax"}, \code{"firstSEmax"} and
  \code{"globalSEmax"}. Default is \code{"firstSEmax"}.
  See \code{'method'} in
  \code{\link[cluster:clusGap]{clusGap}}, Details and
  \emph{Tibshirani et al., 2001} in References.}

  \item{verbose}{logical flag enabling verbose messages to
  track the computing status of the function.  Default is
  \code{TRUE}.}

  \item{clust_trends}{object of class
  \bold{\link{ClusteredTrends}} containing already computed
  trends for the plotted gene sets.  Default is
  \code{NULL}.}

  \item{N_clusters}{an integer that is the number of
  clusters in which the dynamics should be regrouped. The
  cutoff of the clustering tree is automatically calculated
  accordingly.  Default is \code{NULL}, in which case the
  dendrogram is not cut and no clusters are identified.}

  \item{myclusters}{a character vector of colors for
  predefined clusters of the represented genesets, with as
  many levels as the value of \code{N_clusters}. Default is
  \code{NULL}, in which case the clusters are automatically
  identified and colored via the \code{\link{cutree}}
  function and the \code{N_clusters} argument only.}

  \item{label.clusters}{if \code{N_clusters} is not
  \code{NULL}, a character vector of length
  \code{N_clusterss}.  Default is \code{NULL}, in which
  case if \code{N_clusters} is not \code{NULL}, clusters
  are simply labelled with numbers.}

  \item{prev_rowCL}{a \bold{hclust} object, such as the one
  return by the present plotting funstion (see Value) for
  instance.  If not \code{NULL}, no clustering is
  calculated by the present plotting function and this tree
  is used to represent the gene sets dynamics.  Default is
  \code{NULL}.}

  \item{descript}{logical flag indicating that the
  description of the gene sets should appear after their
  name on the right side of the plot if \code{TRUE}.
  Default is \code{TRUE}.  See Details.}

  \item{plotAll}{logical flag indicating wether a first
  heatmap with the median over all the patients should be
  plotted, or not.  Default is \code{TRUE}.}

  \item{color.vec}{a character strings vector used to
  define the color \link[grDevices:palette]{palette} used
  in the plot.  Default is \code{c("#D73027",
  "#FC8D59","lightyellow", "#91BFDB", "#4575B4")}.}

  \item{legend.breaks}{a numeric vector indicating the
  splitting points for coloring.  Default is \code{NULL},
  in which case the break points will be spaced equally and
  symetrically about 0.}

  \item{label.column}{a vector of character strings with
  the labels to be displayed for the columns (i.e. the time
  points).  Default is \code{NULL}.}

  \item{time_unit}{the time unit to be displayed (such as
  \code{"Y"}, \code{"M"}, \code{"W"}, \code{"D"},
  \code{"H"}, etc) next to the values of \code{TimePoint}
  in the columns labels when \code{label.column} is
  \code{NULL}.  Default is \code{""}.}

  \item{cex.label.row}{a numerical value giving the amount
  by which row labels text should be magnified relative to
  the default \code{1}.}

  \item{cex.label.column}{a numerical value giving the
  amount by which column labels text should be magnified
  relative to the default \code{1}.}

  \item{margins}{numeric vector of length 2 containing the
  margins (see \code{\link{par}(mar= *))} for column and
  row names, respectively.  Default is \code{c(15, 100)}.
  See Details.}

  \item{heatKey.size}{the size of the color key for the
  heatmap fill.  Default is \code{1}.}

  \item{dendrogram.size}{the horizontal size of the
  dendrogram.  Default is \code{1}}

  \item{heatmap.height}{the height of the heatmap.  Default
  is \code{1}}

  \item{heatmap.width}{the width of the heatmap.  Default
  is \code{1}}

  \item{cex.clusterKey}{a numerical value giving the amount
  by which the clusters legend text should be magnified
  relative to the default \code{1}, when \code{N_clusters}
  is not \code{NULL}.}

  \item{cex.main}{a numerical value giving the amount by
  which title text should be magnified relative to the
  default \code{1}.}

  \item{horiz.clusterKey}{a logical flag; if \code{TRUE},
  set the legend for clusters horizontally rather than
  vertically.  Only used if the \code{N_clusters} argument
  is not \code{NULL}.  Default is \code{TRUE}.}

  \item{main}{a character string for an optionnal title.
  Default is \code{NULL}.}

  \item{subtitle}{a character string for an optionnal
  subtitle.  Default is \code{NULL}.}

  \item{\dots}{other parameters to be passed through to
  plotting functions.}
}
\value{
  An object of class \bold{\link{hclust}} which describes
  the tree produced by the clustering process.  The object
  is a list with components: \itemize{ \item merge an
  \eqn{n-1} by \eqn{2} matrix.  Row \eqn{i} of \code{merge}
  describes the merging of clusters at step i of the
  clustering. If an element \eqn{j} in the row is negative,
  then observation -\eqn{j} was merged at this stage.  If
  \eqn{j} is positive then the merge was with the cluster
  formed at the (earlier) stage \eqn{j} of the algorithm.
  Thus negative entries in merge indicate agglomerations of
  singletons, and positive entries indicate agglomerations
  of non-singletons. \item height a set of \eqn{n-1} real
  values (non-decreasing for ultrametric trees).  The
  clustering height: that is, the value of the criterion
  associated with the Ward clustering method. \item order a
  vector giving the permutation of the original
  observations suitable for plotting, in the sense that a
  cluster plot using this ordering and matrix merge will
  not have crossings of the branches. \item labels the gene
  sets name. \item call the call which produced the result
  clustering: \cr\code{hclust(d = dist(map2heat, method =
  "euclidean"), method = "ward")} \item method "ward", as
  it is the cluster method that is used. \item dist.method
  "euclidean", as it is the distance that has been used.
  \item legend.breaks a numeric vector giving the splitting
  points used for coloring the heatmap.  If \code{plot} is
  \code{FALSE}, then it is \code{NULL}. \item myclusters a
  character vector of colors for clusters of the
  represented genesets, with as many levels as the value of
  \code{N_clusters}. If no clusters were represented, than
  this is \code{NULL}. \item ddr a \bold{dendrogram} object
  with the reordering used for the heatmap.  See
  \code{\link{heatmap.2}}. \item clustersExport a data
  frame with 2 variables containing the two following
  variables : \itemize{ \item \code{GeneSet}: the gene sets
  clustered.  \item \code{Cluster}: the cluster they belong
  to.  } The data frame is order by the variable
  \code{Cluster}. }
}
\description{
  This function plots a series of gene sets dynamic trends
  heatmaps.  One heatmap is drawned for each patient.
}
\details{
  On the heatmap, each line corresponds to a gene set, and
  each column to a timepoint.

  First a heatmap is computed on all the patients (see
  \code{\link{plot.TcGSA}} and \code{\link{clustTrend}}) to
  define the clustering. Then, the clustering and coloring
  thus defined on all the patients are consistently used in
  the separate heatmaps that are plotted by patient.

  If \code{expr} is a matrix or a dataframe, then the
  "original" data are plotted.  On the other hand, if
  \code{expr} is a list returned in the
  \code{'Estimations'} element of \code{\link{TcGSA.LR}},
  then it is those "estimations" made by the
  \code{\link{TcGSA.LR}} function that are plotted.

  If \code{descript} is \code{FALSE}, the second element of
  \code{margins} can be reduced (for instance use
  \code{margins = c(5, 10)}), as there is not so much need
  for space in order to display only the gene set names,
  without their description.

  The median shown in the heatmap uses the respectively
  standardized (reduced and centered) expression of the
  genes over the patients.
}
\examples{
data(data_simu_TcGSA)

tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design,
													 subject_name="Patient_ID", time_name="TimePoint",
                          time_func="linear", crossedRandom=FALSE)

plotPat.TcGSA(x=tcgsa_sim_1grp, expr=expr_1grp,
    Subject_ID=design$Patient_ID, TimePoint=design$TimePoint,
    B=100,
    time_unit="H"
    )

plotPat.TcGSA(x=tcgsa_sim_1grp, expr=tcgsa_sim_1grp$Estimations,
    Subject_ID=design$Patient_ID, TimePoint=design$TimePoint,
    baseline=1,
    B=100,
    time_unit="H"
    )
}
\author{
  Boris P. Hejblum
}
\references{
  Hejblum, B.P., Skinner, J., Thiebaut, R., 2013, TcGSA: a
  gene set approach for longitudinal gene expression data
  analysis, \bold{submitted}.
}
\seealso{
  \code{\link{plot.TcGSA}},
  \code{\link[gplots:heatmap.2]{heatmap.2}},
  \code{\link{TcGSA.LR}}, \code{\link{hclust}}
}

