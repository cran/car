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

\method{scatterplotMatrix}{formula}(x, data=NULL, subset, labels, ...)

\method{scatterplotMatrix}{default}(x, var.labels = colnames(x), diagonal = c("density", 
    "boxplot", "histogram", "oned", "qqplot", "none"), adjust = 1, 
    nclass, plot.points = TRUE, smooth = TRUE, 
    spread = smooth && !by.groups, span = 0.5, 
    loess.threshold = 5, reg.line = lm, 
    transform = FALSE, family = c("bcPower", "yjPower"), ellipse = FALSE, 
    levels = c(0.5, 0.95), robust = TRUE, groups = NULL, by.groups = FALSE, 
    labels, id.method="mahal", id.n=0, id.cex=1, id.col=palette()[1], 
    col = if (n.groups == 1) palette()[2:1] else rep(palette(), 
        length = n.groups), 
    pch = 1:n.groups, lwd = 1, lwd.smooth = lwd, 
    lwd.spread = lwd, lty = 1, lty.smooth = lty, lty.spread = 2, 
    cex = par("cex"), cex.axis = par("cex.axis"), cex.labels = NULL, 
    cex.main = par("cex.main"), legend.plot = length(levels(groups)) > 
        1, row1attop = TRUE, ...)

spm(x, ...)
}

\arguments{
  \item{x}{a data matrix, numeric data frame, or a one-sided \dQuote{model} formula, of the form 
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
  \item{diagonal}{contents of the diagonal panels of the plot.}
  \item{adjust}{relative bandwidth for density estimate, passed to 
    \code{density} function.}
  \item{nclass}{number of bins for histogram, passed to \code{hist}
    function.}
  \item{plot.points}{if \code{TRUE} the points are plotted in each 
    off-diagonal panel.}
  \item{smooth}{if \code{TRUE} a loess smooth is plotted in each
    off-diagonal panel.}
  \item{spread}{if \code{TRUE} (the default when not smoothing by groups), a smoother is applied to
  	the root-mean-square positive and negative residuals from the loess line to display conditional spread and
  	asymmetry.}
  \item{span}{span for loess smoother.}
   \item{loess.threshold}{suppress the loess smoother if there are fewer than \code{loess.threshold} unique
  	values (default, 5) of the variable on the vertical axis.}
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
  \item{pch}{plotting characters for points; default is the plotting characters in
    order (see \code{\link[graphics]{par}}).}
  \item{col}{colors for lines and points; the default is taken from the color palette, 
  	with \code{palette()[2]} for nonparametric regression lines and \code{palette()[1]} for
  	linear regression lines and points if there are
  	no groups, and successive colors for the groups if there are groups.}
  \item{lwd}{width of linear-regression lines (default \code{1}).}
  \item{lwd.smooth}{width for smooth regression lines (default is the same as \code{lwd}).}
  \item{lwd.spread}{width for lines showing spread (default is the same as \code{lwd}).}
  \item{lty}{type of linear-regression lines (default \code{1}, solid line).}
  \item{lty.smooth}{type of smooth regression lines (default is the same as \code{lty}).}
  \item{lty.spread}{width for lines showing spread (default is \code{2}, broken line).}
  \item{cex, cex.axis, cex.labels, cex.main}{set sizes of various graphical elements
    (see \code{\link[graphics]{par}}).}
  \item{legend.plot}{if \code{TRUE} then a legend for the groups is plotted
    in the first diagonal cell.}
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
  \code{\link{showLabels}}.}

\examples{
scatterplotMatrix(~ income + education + prestige | type, data=Duncan)
scatterplotMatrix(~ income + education + prestige, 
    transform=TRUE, data=Duncan)
scatterplotMatrix(~ income + education + prestige | type, smooth=FALSE, 
	by.group=TRUE, transform=TRUE, data=Duncan)
}

\keyword{hplot}
