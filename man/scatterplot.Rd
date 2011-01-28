\name{scatterplot}
\alias{scatterplot}
\alias{scatterplot.formula}
\alias{scatterplot.default}
\alias{sp}

\title{Scatterplots with Boxplots}
\description{
  Makes enhanced scatterplots, with boxplots in the margins, a lowess
  smooth, smoothed conditional spread, outlier identification, and a regression line; \code{sp} is an abbreviation for
  \code{scatterplot}.
}
\usage{
scatterplot(x, ...)

\method{scatterplot}{formula}(x, data, subset, xlab, ylab, legend.title, legend.coords,
	labels, ...)

\method{scatterplot}{default}(x, y, smooth = TRUE, 
    spread = !by.groups, span = 0.5, 
    loess.threshold = 5, reg.line = lm, 
    boxplots = if (by.groups) "" else "xy", 
    xlab = deparse(substitute(x)), ylab = deparse(substitute(y)), 
    las = par("las"), lwd = 1, lwd.smooth = lwd, lwd.spread = lwd, 
    lty = 1, lty.smooth = lty, lty.spread = 2, labels, 
    id.method = "mahal", 
    id.n = if(id.method[1]=="identify") length(x) else 0, 
    id.cex = 1, id.col = palette()[1], log = "", jitter = list(), 
    xlim = NULL, ylim = NULL, cex = par("cex"), 
    cex.axis = par("cex.axis"), 
    cex.lab = par("cex.lab"), cex.main = par("cex.main"), 
    cex.sub = par("cex.sub"), groups, by.groups = !missing(groups), 
    legend.title = deparse(substitute(groups)), legend.coords,
    ellipse = FALSE, levels = c(0.5, 0.95), robust = TRUE, 
    col = if (n.groups == 1) palette()[2:1] else rep(palette(), 
    length = n.groups), 
    pch = 1:n.groups, legend.plot = !missing(groups), reset.par = TRUE, 
    grid=TRUE, ...)
  
sp(...)
}

\arguments{
  \item{x}{vector of horizontal coordinates, or a ``model'' formula, of the form \code{y ~ x} or 
    (to plot by groups) \code{y ~ x | z}, where \code{z} evaluates to a factor 
    or other variable dividing the data into groups. If \code{x} is a factor, then parallel boxplots
    are produced using the \code{\link{Boxplot}} function.}
  \item{y}{vector of vertical coordinates.}
  \item{data}{data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{smooth}{if \code{TRUE} (the default) a loess nonparametric regression line is
    drawn on the plot.}
  \item{spread}{if \code{TRUE} (the default when there are no groups), a smoother is applied to
  	the root-mean-square positive and negative residuals from the loess line to display conditional spread and
  	asymmetry.}
  \item{span}{span for the loess smoother.}
  \item{loess.threshold}{suppress the loess smoother if there are fewer than \code{loess.threshold} unique
  	values (default, 5) of \code{y}.}
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
  \item{lwd.smooth}{width for smooth regression lines (default is the same as \code{lwd}).}
  \item{lwd.spread}{width for lines showing spread (default is the same as \code{lwd}).}
  \item{lty}{type of linear-regression lines (default \code{1}, solid line).}
  \item{lty.smooth}{type of smooth regression lines (default is the same as \code{lty}).}
  \item{lty.spread}{width for lines showing spread (default is \code{2}, broken line).}
  \item{id.method,id.n,id.cex,id.col}{Arguments for the labelling of
    points.  The default is \code{id.n=0} for labeling no points.  See
    \code{\link{showLabels}} for details of these arguments.  If the plot uses
    different colors for groups, then the \code{id.col} argument is ignored and
    label colors are determined by the \code{col} argument. 
}
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
  	with \code{palette()[2]} for nonparametric regression lines and \code{palette()[1]} 
  	for linear regression line and points if there are
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
  \code{\link{showLabels}}.}

\examples{
scatterplot(prestige ~ income, data=Prestige, ellipse=TRUE)

scatterplot(prestige ~ income|type, data=Prestige, span=1)

scatterplot(prestige ~ income|type, data=Prestige, span=1, 
	legend.coords="topleft")

scatterplot(vocabulary ~ education, jitter=list(x=1, y=1), 
	data=Vocab, id.n=0, smooth=FALSE)

scatterplot(infant.mortality ~ gdp, log="xy", data=UN, id.n=5)

scatterplot(income ~ type, data=Prestige)

	\dontrun{
scatterplot(infant.mortality ~ gdp, id.method="identify", data=UN)
	}
}

\keyword{hplot}
