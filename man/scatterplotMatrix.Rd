\name{scatterplotMatrix}
\alias{scatterplotMatrix}
\alias{scatterplotMatrix.formula}
\alias{scatterplotMatrix.default}
\alias{spm}
\title{Scatterplot Matrices}
\description{
  Enhanced scatterplot matrices with univariate displays down the diagonal;
  \code{spm} is an abbreviation for \code{scatterplotMatrix}.
  This function just sets up a call to \code{pairs} with custom panel functions.
}
\usage{

scatterplotMatrix(x, ...)

\method{scatterplotMatrix}{formula}(formula, data=NULL, subset, labels, ...)

\method{scatterplotMatrix}{default}(x, var.labels=colnames(x), 
    diagonal=c("density", "boxplot", "histogram", "oned", "qqplot", "none"), 
    adjust=1, nclass,
    plot.points=TRUE, smoother=loessLine, smoother.args=list(), smooth, span,
    spread = !by.groups, reg.line=lm,
    transform=FALSE, family=c("bcPower", "yjPower"),
    ellipse=FALSE, levels=c(.5, .95), robust=TRUE,
    groups=NULL, by.groups=FALSE, 
    use=c("complete.obs", "pairwise.complete.obs"),
    labels, id.method="mahal", id.n=0, id.cex=1, id.col=palette()[1],
    col=if (n.groups == 1) palette()[3:1] else rep(palette(), length=n.groups),
    pch=1:n.groups, lwd=1, lty=1, 
    cex=par("cex"), cex.axis=par("cex.axis"), cex.labels=NULL, 
    cex.main=par("cex.main"), 
    legend.plot=length(levels(groups)) > 1, legend.pos=NULL, row1attop=TRUE, ...)

spm(x, ...)
}

\arguments{
  \item{x}{a data matrix, numeric data frame.}
  \item{formula}{a one-sided \dQuote{model} formula, of the form 
    \code{ ~ x1 + x2 + ... + xk} or \code{ ~ x1 + x2 + ... + xk | z} where \code{z}
    evaluates to a factor or other variable to divide the data into groups.}
  \item{data}{for \code{scatterplotMatrix.formula},
    a data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{labels,id.method,id.n,id.cex,id.col}{Arguments for the labelling of
    points.  The default is \code{id.n=0} for labeling no points.  See
    \code{\link{showLabels}} for details of these arguments.  If the plot uses
    different colors for groups, then the \code{id.col} argument is ignored and
    label colors are determined by the \code{col} argument. 
}
  \item{var.labels}{variable labels (for the diagonal of the plot).}
  \item{diagonal}{contents of the diagonal panels of the plot. If plotting by groups, a different
    univariate display (with the exception of \code{"histogram"}) will be drawn for each group.}
  \item{adjust}{relative bandwidth for density estimate, passed to 
    \code{density} function.}
  \item{nclass}{number of bins for histogram, passed to \code{hist}
    function.}
  \item{plot.points}{if \code{TRUE} the points are plotted in each 
    off-diagonal panel.}
  \item{smoother}{a function to draw a nonparametric-regression smooth; the default is \code{\link{gamLine}}, which
    uses the \code{\link[mgcv]{gam}} function in the \pkg{mgcv} package. For this and other smoothers,
    see \code{\link{ScatterplotSmoothers}}.
    Setting this argument to something other than a function, e.g., \code{FALSE} suppresses the smoother.}
  \item{smoother.args}{a list of named values to be passed to the smoother function; the specified elements of the
    list depend upon the smoother (see \code{\link{ScatterplotSmoothers}}).}
  \item{smooth, span}{these arguments are included for backwards compatility: if \code{smooth=TRUE} then \code{smoother} is set to \code{loessLine},
  and if \code{span} is specified, it is added to \code{smoother.args}.}
  \item{spread}{if TRUE, estimate the (square root) of the variance function.  For \code{loessLine} and
    for \code{gamLine}, this is done by separately smoothing the squares of the postive and negative
    residuals from the mean fit, and then adding the square root of the fitted values to the mean fit.  For
    \code{quantregLine}, fit the .25 and .75 quantiles with a quantile regression additive model.
    The default is TRUE if \code{by.groups=FALSE} and FALSE is \code{by.groups=TRUE}.}
  \item{reg.line}{if not \code{FALSE} a line is plotted using the
    function given by this argument; e.g., using \code{rlm} in
    package \code{MASS} plots a robust-regression line.}
  \item{transform}{if \code{TRUE}, multivariate normalizing power transformations
    are computed with \code{\link{powerTransform}}, rounding the estimated powers to `nice' values for plotting; 
    if a vector of powers, one for each variable, these are applied prior to plotting. If  there are \code{groups}
    and \code{by.groups} is \code{TRUE}, then the transformations are estimated \emph{conditional} on the
    \code{groups} factor.}
  \item{family}{family of transformations to estimate: \code{"bcPower"} for the Box-Cox family or 
    \code{"yjPower"} for the Yeo-Johnson family (see \code{\link{powerTransform}}).}
  \item{ellipse}{if \code{TRUE} data-concentration ellipses are plotted in
    the off-diagonal panels.}
  \item{levels}{levels or levels at which concentration ellipses are plotted;
    the default is \code{c(.5, .9)}.}
  \item{robust}{if \code{TRUE} use the \code{cov.trob} function in the \code{MASS} package
    to calculate the center and covariance matrix for the data ellipses.}
  \item{groups}{a factor or other variable dividing the data into groups; groups are
    plotted with different colors and plotting characters.}
  \item{by.groups}{if \code{TRUE}, regression lines are fit by groups.}
  \item{use}{if \code{"complete.obs"} (the default), cases with missing data are omitted; if \code{"pairwise.complete.obs"), all valid cases are used
    in each panel of the plot.}}
  \item{pch}{plotting characters for points; default is the plotting characters in
    order (see \code{\link[graphics]{par}}).}
  \item{col}{colors for lines and points; the default is taken from the color palette, 
  	with \code{palette()[3]} for linear regression lines, \code{palette()[2]} for nonparametric regression lines, and \code{palette()[1]} 
  	for points if there are
  	no groups, and successive colors for the groups if there are groups.}
  \item{lwd}{width of linear-regression lines (default \code{1}).}
  \item{lty}{type of linear-regression lines (default \code{1}, solid line).}
  \item{cex, cex.axis, cex.labels, cex.main}{set sizes of various graphical elements
    (see \code{\link[graphics]{par}}).}
  \item{legend.plot}{if \code{TRUE} then a legend for the groups is plotted
    in the first diagonal cell.}
  \item{legend.pos}{position for the legend, specified as one of the keywords accepted by
    \code{\link{legend}}. If \code{NULL}, the default, the position will vary by the 
    \code{diagonal} argument --- e.g., \code{"topright"} for \code{diagonal="density"}.}
  \item{row1attop}{If \code{TRUE} (the default) the first row is at the top, as in a matrix, as
  	opposed to at the bottom, as in graph (argument suggested by Richard Heiberger).}
  \item{...}{arguments to pass down.}
}

\value{
  \code{NULL}. This function is used for its side effect: producing
  a plot.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}

\seealso{\code{\link[graphics]{pairs}}, \code{\link{scatterplot}},
  \code{\link{dataEllipse}}, \code{\link{powerTransform}}, 
  \code{\link{bcPower}}, \code{\link{yjPower}}, \code{\link[MASS]{cov.trob}},
  \code{\link{showLabels}}, \code{\link{ScatterplotSmoothers}}.}

\examples{
scatterplotMatrix(~ income + education + prestige | type, data=Duncan)
scatterplotMatrix(~ income + education + prestige, 
    transform=TRUE, data=Duncan, smoother=loessLine)
scatterplotMatrix(~ income + education + prestige | type, smoother=FALSE, 
	by.group=TRUE, transform=TRUE, data=Duncan)
}

\keyword{hplot}
