\name{scatterplot}
\alias{scatterplot}
\alias{scatterplot.formula}
\alias{scatterplot.default}
\alias{sp}

\title{Scatterplots with Boxplots}
\description{
  Makes enhanced scatterplots, with boxplots in the margins, a nonparametric regression smooth, 
  smoothed conditional spread, outlier identification, and a regression line; 
  \code{sp} is an abbreviation for \code{scatterplot}.
}
\usage{
scatterplot(x, ...)

\method{scatterplot}{formula}(formula, data, subset, xlab, ylab, legend.title, legend.coords,
	labels, ...)

\method{scatterplot}{default}(x, y, 
    smoother=loessLine, smoother.args=list(), smooth, span,
    spread=!by.groups, reg.line=lm,
    boxplots=if (by.groups) "" else "xy",
    xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), las=par("las"),
    lwd=1, lty=1,
    labels, id.method = "mahal", 
    id.n = if(id.method[1]=="identify") length(x) else 0, 
    id.cex = 1, id.col = palette()[1],
    log="", jitter=list(), xlim=NULL, ylim=NULL,
    cex=par("cex"), cex.axis=par("cex.axis"), cex.lab=par("cex.lab"), 
    cex.main=par("cex.main"), cex.sub=par("cex.sub"), 
    groups, by.groups=!missing(groups), 
    legend.title=deparse(substitute(groups)), legend.coords,
    ellipse=FALSE, levels=c(.5, .95), robust=TRUE,
    col=if (n.groups == 1) palette()[3:1] else rep(palette(), length=n.groups),
    pch=1:n.groups, 
    legend.plot=!missing(groups), reset.par=TRUE, grid=TRUE, ...)

sp(...)
}

\arguments{
  \item{x}{vector of horizontal coordinates}
  \item{y}{vector of vertical coordinates.}
  \item{formula}{a ``model'' formula, of the form \code{y ~ x} or 
    (to plot by groups) \code{y ~ x | z}, where \code{z} evaluates to a factor 
    or other variable dividing the data into groups. If \code{x} is a factor, then parallel boxplots
    are produced using the \code{\link{Boxplot}} function.}
  \item{data}{data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{smoother}{a function to draw a nonparametric-regression smooth; the default is \code{\link{loessLine}},
    which does loess smoothing.  The function \code{\link{gamLine}} fits a generalized additive model and
    allows including a link and error function.
    See \code{\link{ScatterplotSmoothers}}.
    Setting this argument to something other than a function, e.g., \code{FALSE} suppresses the smoother.}
  \item{smoother.args}{a list of named values to be passed to the smoother function; the specified elements of the
    list depend upon the smoother (see \code{\link{ScatterplotSmoothers}}).}
  \item{smooth, span}{these arguments are included for backwards compatility: 
    if \code{smooth=TRUE} then \code{smoother} is set to \code{loessLine},
    and if \code{span} is specified, it is added to \code{smoother.args}.}
  \item{spread}{if TRUE, estimate the (square root) of the variance function.  For \code{loessLine} and
    for \code{gamLine}, this is done by separately smoothing the squares of the postive and negative
    residuals from the mean fit, and then adding the square root of the fitted values to the mean fit.  For
    \code{quantregLine}, fit the .25 and .75 quantiles with a quantile regression additive model.
    The default is TRUE if \code{by.groups=FALSE} and FALSE is \code{by.groups=TRUE}.}
  \item{reg.line}{function to draw a regression line on the plot
    or \code{FALSE} not to plot a regression line.}
  \item{boxplots}{if \code{"x"} a boxplot for \code{x} is drawn below the plot;
    if \code{"y"} a boxplot for \code{y} is drawn to the left of the plot;
    if \code{"xy"} both boxplots are drawn; set to \code{""} or \code{FALSE} to
    suppress both boxplots.}
  \item{xlab}{label for horizontal axis.}
  \item{ylab}{label for vertical axis.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{lwd}{width of linear-regression lines (default \code{1}).}
  \item{lty}{type of linear-regression lines (default \code{1}, solid line).}
  \item{id.method,id.n,id.cex,id.col}{Arguments for the labelling of
    points.  The default is \code{id.n=0} for labeling no points.  See
    \code{\link{showLabels}} for details of these arguments.  If the plot uses
    different colors for groups, then the \code{id.col} argument is ignored and
    label colors are determined by the \code{col} argument.}
  \item{labels}{a vector of point labels; if absent, the function tries to determine reasonable
  	labels, and, failing that, will use observation numbers.}
  \item{log}{same as the \code{log} argument to \code{\link{plot}}, to produce log axes.}
  \item{jitter}{a list with elements \code{x} or \code{y} or both, specifying jitter factors
    for the horizontal and vertical coordinates of the points in the scatterplot. The 
    \code{\link[base]{jitter}} function is used to randomly perturb the points; specifying a
    factor of \code{1} produces the default jitter. Fitted lines are unaffected by the jitter.}
  \item{xlim}{the x limits (min, max) of the plot; if \code{NULL}, determined from the data.}
  \item{ylim}{the y limits (min, max) of the plot; if \code{NULL}, determined from the data.}
  \item{groups}{a factor or other variable dividing the data into groups; groups are
    plotted with different colors and plotting characters.}
  \item{by.groups}{if \code{TRUE}, regression lines are fit by groups.}
  \item{legend.title}{title for legend box; defaults to the name of the groups variable.}
  \item{legend.coords}{coordinates for placing legend; an be a list with components x and y to
  	specify the coordinates of the upper-left-hand corner of the legend; or a quoted keyword, such as \code{"topleft"},
  	recognized by \code{\link{legend}}.}
  \item{ellipse}{if \code{TRUE} data-concentration ellipses are plotted.}
  \item{levels}{level or levels at which concentration ellipses are plotted;
    the default is \code{c(.5, .95)}.}
  \item{robust}{if \code{TRUE} (the default) use the \code{cov.trob} function in the \code{MASS} package
    to calculate the center and covariance matrix for the data ellipses.}
  \item{col}{colors for lines and points; the default is taken from the color palette, 
  	with \code{palette()[3]} for linear regression lines, \code{palette()[2]} for nonparametric regression lines, and \code{palette()[1]} 
  	for points if there are
  	no groups, and successive colors for the groups if there are groups.}
  \item{pch}{plotting characters for points; default is the plotting characters in
    order (see \code{\link[graphics]{par}}).}
  \item{cex, cex.axis, cex.lab, cex.main, cex.sub}{set sizes of various graphical elements; 
    (see \code{\link[graphics]{par}}).}
  \item{legend.plot}{if \code{TRUE} then a legend for the groups is plotted in the upper margin.}
  \item{reset.par}{if \code{TRUE} then plotting parameters are reset to their previous values
    when \code{scatterplot} exits; if \code{FALSE} then the \code{mar} and \code{mfcol} parameters are
    altered for the current plotting device. Set to \code{FALSE} if you want to add graphical elements
    (such as lines) to the plot.}
  \item{\dots}{other arguments passed down and to \code{plot}.}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
}

\value{
  If points are identified, their labels are returned; otherwise \code{NULL} is returned invisibly.
}

\author{John Fox \email{jfox@mcmaster.ca}}


\seealso{\code{\link[graphics]{boxplot}},
  \code{\link[base]{jitter}}, \code{\link{legend}},
  \code{\link{scatterplotMatrix}}, \code{\link{dataEllipse}}, \code{\link{Boxplot}},
  \code{\link[MASS]{cov.trob}},
  \code{\link{showLabels}}, \code{\link{ScatterplotSmoothers}}.}

\examples{
scatterplot(prestige ~ income, data=Prestige, ellipse=TRUE)

if (interactive()){
scatterplot(prestige ~ income, data=Prestige, smoother=quantregLine)
}

scatterplot(prestige ~ income|type, data=Prestige, smoother=loessLine, 
    smoother.args=list(span=1))

scatterplot(prestige ~ income|type, data=Prestige, legend.coords="topleft")

scatterplot(vocabulary ~ education, jitter=list(x=1, y=1), 
	data=Vocab, id.n=0, smoother=FALSE)

scatterplot(infant.mortality ~ gdp, log="xy", data=UN, id.n=5)

scatterplot(income ~ type, data=Prestige)

	\dontrun{
scatterplot(infant.mortality ~ gdp, id.method="identify", data=UN)

scatterplot(infant.mortality ~ gdp, id.method="identify", smoother=loessLine, data=UN)
	}
}

\keyword{hplot}
