\name{plot1GS}
\alias{plot1GS}
\title{Plotting a Specific Gene Set}
\usage{
  plot1GS(expr, gmt, Subject_ID, TimePoint, geneset.name,
    baseline = NULL, group.var = NULL,
    Group_ID_paired = NULL, ref = NULL,
    group_of_interest = NULL, FUNcluster = NULL,
    clustering_metric = "euclidian",
    clustering_method = "ward", B = 500, max_trends = 4,
    aggreg.fun = "median", trend.fun = "median",
    methodOptiClust = "firstSEmax", indiv = "genes",
    verbose = TRUE, clustering = TRUE, showTrend = TRUE,
    smooth = TRUE, time_unit = "", title = NULL,
    y.lab = NULL, desc = TRUE, lab.cex = 1, axis.cex = 1,
    main.cex = 1, y.lab.angle = 90, x.axis.angle = 45,
    y.lim = NULL, x.lim = NULL, gg.add = list(theme()),
    plot = TRUE)
}
\arguments{
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

  \item{gmt}{a \bold{gmt} object containing the gene sets
  definition.  See
  \code{\link[GSA:GSA.read.gmt]{GSA.read.gmt}} and
  definition on
  \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29}{www.broadinstitute.org}.}

  \item{Subject_ID}{a factor of length \eqn{p} that is in
  the same order as the columns of \code{expr} (when it is
  a dataframe) and that contains the patient identifier of
  each sample.}

  \item{TimePoint}{a numeric vector or a factor of length
  \eqn{p} that is in the same order as \code{TimePoint} and
  the columns of \code{expr} (when it is a dataframe), and
  that contains the time points at which gene expression
  was measured.}

  \item{geneset.name}{a character string containing the
  name of the gene set to be plotted, that must appear in
  the \code{"geneset.names"} element of \code{gmt}.}

  \item{baseline}{a character string which is the value of
  \code{TimePoint} that can be used as a baseline.  Default
  is \code{NULL}, in which case no timepoint is used as a
  baseline value for gene expression.  Has to be
  \code{NULL} when comparing two treatment groups.  See
  Details.}

  \item{group.var}{in the case of several treatment groups,
  this is a factor of length \eqn{p} that is in the same
  order as \code{Timepoint}, \code{Subject_ID} and the
  columns of \code{expr}.  It indicates to which treatment
  group each sample belongs to.  Default is \code{NULL},
  which means that there is only one treatment group.  See
  Details.}

  \item{Group_ID_paired}{a character vector of length
  \eqn{p} that is in the same order as \code{Timepoint},
  \code{Subject_ID}, \code{group.var} and the columns of
  \code{expr}.  This argument must not be \code{NULL} in
  the case of a paired analysis, and must be \code{NULL}
  otherwise.  Default is \code{NULL}.  See Details.}

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

  \item{FUNcluster}{a function which accepts as first
  argument a matrix \code{x} and as second argument the
  number of clusters desired \code{k}, and which returns a
  list with a component named \code{'cluster'} which is a
  vector of length n = nrow(x) of integers in 1:k,
  determining the clustering or grouping of the n
  observations.  Default is \code{NULL}, in which case a
  hierachical clustering is performed via the function
  \code{\link[cluster:agnes]{agnes}}, using the metric
  \code{clustering_metric} and the method
  \code{clustering_method}.  See \code{'FUNcluster'} in
  \code{\link[cluster:clusGap]{clusGap}} and Details.}

  \item{clustering_metric}{character string specifying the
  metric to be used for calculating dissimilarities between
  observations in the hierarchical clustering when
  \code{FUNcluster} is \code{NULL}.  The currently
  available options are \code{"euclidean"} and
  \code{"manhattan"}.  Default is \code{"euclidean"}.  See
  \code{\link[cluster:agnes]{agnes}}.  Also, a \code{"sts"}
  option is available in TcGSA.  It implements the 'Short
  Time Series' distance [Möller-Levet et al., Fuzzy
  CLustering of short time series and unevenly distributed
  sampling points, \emph{Advances in Intelligent Data
  Analysis V}:330-340 Springer, 2003] designed specifically
  for clustering time series.}

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
  defined statistics function that returns a single numeric
  value.  It specifies the function used to aggregate the
  observations before the clustering.  Default is to
  \code{median}.}

  \item{trend.fun}{a character string such as
  \code{"mean"}, \code{"median"} or the name of any other
  function that returns a single numeric value.  It
  specifies the function used to calculate the trends of
  the identified clustered.  Default is to \code{median}.}

  \item{methodOptiClust}{character string indicating how
  the "optimal" number of clusters is computed from the gap
  statistics and their standard deviations. Possible values
  are \code{"globalmax"}, \code{"firstmax"},
  \code{"Tibs2001SEmax"}, \code{"firstSEmax"} and
  \code{"globalSEmax"}. Default is \code{"firstSEmax"}.
  See \code{'method'} in
  \code{\link[cluster:clusGap]{clusGap}}, Details and
  \emph{Tibshirani et al., 2001} in References.}

  \item{indiv}{a character string indicating by which unit
  observations are aggregated (through \code{aggreg.fun})
  before the clustering.  Possible values are
  \code{"genes"} or \code{"patients"}.  Default is
  \code{"genes"}. See Details.}

  \item{verbose}{logical flag enabling verbose messages to
  track the computing status of the function.  Default is
  \code{TRUE}.}

  \item{clustering}{logical flag.  If \code{FALSE}, there
  is no clustering representation; if \code{TRUE}, the
  lines are colored according to which cluster they belong
  to.  Default is \code{TRUE}.  See Details.}

  \item{showTrend}{logical flag.  If \code{TRUE}, a black
  line is added for each cluster, representing the
  corresponding \code{trend.fun}.  Default is \code{TRUE}.}

  \item{smooth}{logical flag.  If \code{TRUE} and
  \code{showTrend} is also \code{TRUE}, the representation
  of each cluster \code{trend.fun} is smoothed using cubic
  polynoms (see
  \code{\link[ggplot2:stat_smooth]{stat_smooth}}. Default
  is \code{TRUE}.}

  \item{time_unit}{the time unit to be displayed (such as
  \code{"Y"}, \code{"M"}, \code{"W"}, \code{"D"},
  \code{"H"}, etc) next to the values of \code{TimePoint}
  on the x-axis.  Default is \code{""}, in which case the
  time scale on the x-axis is proportionnal to the time
  values.}

  \item{title}{character specifying the title of the plot.
  If \code{NULL}, a title is automatically generated, if
  \code{""}, no title appears.  Default is \code{NULL}.}

  \item{y.lab}{character specifying the annotation of the y
  axis.  If \code{NULL}, an annotation is automatically
  generated, if \code{""}, no annotation appears.  Default
  is \code{NULL}.}

  \item{desc}{a logical flag. If \code{TRUE}, a line is
  added to the title of the plot with the description of
  the gene set plotted (from the gmt file). Default is
  \code{TRUE}.}

  \item{lab.cex}{a numerical value giving the amount by
  which lab labels text should be magnified relative to the
  default \code{1}.}

  \item{axis.cex}{a numerical value giving the amount by
  which axis annotation text should be magnified relative
  to the default \code{1}.}

  \item{main.cex}{a numerical value giving the amount by
  which title text should be magnified relative to the
  default \code{1}.}

  \item{y.lab.angle}{a numerical value (in [0, 360]) giving
  the orientation by which y-label text should be turned
  (anti-clockwise).  Default is \code{90}. See
  \code{\link[ggplot2:element_text]{element_text}}.}

  \item{x.axis.angle}{a numerical value (in [0, 360])
  giving the orientation by which x-axis annotation text
  should be turned (anti-clockwise).  Default is
  \code{45}.}

  \item{y.lim}{a numeric vector of length 2 giving the
  range of the y-axis. See \code{\link{plot.default}}.}

  \item{x.lim}{if numeric, will create a continuous scale,
  if factor or character, will create a discrete scale.
  Observations not in this range will be dropped.  See
  \code{\link{xlim}}.}

  \item{gg.add}{A list of instructions to add to the
  ggplot2 instruction.  See \link{+.gg}.  Default is
  \code{list(theme())}, which adds nothing to the plot.}

  \item{plot}{logical flag.  If \code{FALSE}, no plot is
  drawn.  Default is \code{TRUE}.}
}
\value{
  A dataframe the 2 following variables: \itemize{ \item
  \code{ProbeID} which contains the IDs of the probes of
  the plotted gene set. \item \code{Cluster} which to which
  cluster the probe belongs to. } If \code{clustering} is
  \code{FALSE}, then \code{Cluster} is \code{NA} for all
  the probes.
}
\description{
  This function can plot different representations of the
  gene expression in a specific gene set.
}
\details{
  If \code{expr} is a matrix or a dataframe, then the
  "original" data are plotted.  On the other hand, if
  \code{expr} is a list returned in the
  \code{'Estimations'} element of \code{\link{TcGSA.LR}},
  then it is those "estimations" made by the
  \code{\link{TcGSA.LR}} function that are plotted.

  If \code{indiv} is 'genes', then each line of the plot is
  the median of a gene expression over the patients. On the
  other hand, if \code{indiv} is 'patients', then each line
  of the plot is the median of a patient genes expression
  in this gene set.

  This function uses the Gap statistics to determine the
  optimal number of clusters in the plotted gene set.  See
  \code{\link[cluster:clusGap]{clusGap}}.
}
\examples{
data(data_simu_TcGSA)
tcgsa_sim_1grp <- TcGSA.LR(expr=expr_1grp, gmt=gmt_sim, design=design,
													 subject_name="Patient_ID", time_name="TimePoint",
                          time_func="linear", crossedRandom=FALSE)

plot1GS(expr=expr_1grp, TimePoint=design$TimePoint, Subject_ID=design$Patient_ID, gmt=gmt_sim,
       geneset.name="Gene set 4",
       indiv="genes", clustering=FALSE,
       time_unit="H",
       lab.cex=0.7)

plot1GS(expr=expr_1grp, TimePoint=design$TimePoint, Subject_ID=design$Patient_ID, gmt=gmt_sim,
       geneset.name="Gene set 4",
       indiv="patients", clustering=FALSE, baseline=1,
       time_unit="H",
       lab.cex=0.7)

plot1GS(expr=tcgsa_sim_1grp$Estimations, TimePoint=design$TimePoint, Subject_ID=design$Patient_ID, gmt=gmt_sim,
       geneset.name="Gene set 3",
       indiv="genes",
       time_unit="H",
       lab.cex=0.7
)


colval <- c(hsv(0.56, 0.9, 1),
           hsv(0, 0.27, 1),
           hsv(0.52, 1, 0.5),
           hsv(0, 0.55, 0.97),
           hsv(0.66, 0.15, 1),
           hsv(0, 0.81, 0.55),
           hsv(0.7, 1, 0.7),
           hsv(0.42, 0.33, 1)
)
n <- length(colval);  y <- 1:n
op <- par(mar=rep(1.5,4))
plot(y, axes = FALSE, frame.plot = TRUE,
		 xlab = "", ylab = "", pch = 21, cex = 8,
		 bg = colval, ylim=c(-1,n+1), xlim=c(-1,n+1),
		 main = "Color scale"
)
par(op)


plot1GS(expr=expr_1grp, TimePoint=design$TimePoint, Subject_ID=design$Patient_ID, gmt=gmt_sim,
       geneset.name="Gene set 5",
       indiv="genes",
       time_unit="H",
       title="",
       gg.add=list(scale_color_manual(values=colval), guides(colour = guide_legend(reverse=TRUE))),
       lab.cex=0.7
)
}
\author{
  Boris P. Hejblum
}
\references{
  Tibshirani, R., Walther, G. and Hastie, T., 2001,
  Estimating the number of data clusters via the Gap
  statistic, \emph{Journal of the Royal Statistical
  Society, Series B (Statistical Methodology)}, \bold{63},
  2: 411--423.
}
\seealso{
  \code{\link{ggplot2}},
  \code{\link[cluster:clusGap]{clusGap}}
}

