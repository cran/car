\name{scatterplot.matrix}
\alias{scatterplot.matrix}
\alias{scatterplot.matrix.formula}
\alias{scatterplot.matrix.default}
\alias{spm}
\title{Scatterplot Matrices}
\description{
  Scatterplot matrices with univariate displays down the diagonal;
  \code{spm} is an abbreviation for \code{scatterplot.matrix}.
  This function just sets up a call to \code{pairs}.
}
\usage{

scatterplot.matrix(x, ...)

\method{scatterplot.matrix}{formula}(formula, data=NULL, subset, ...)

\method{scatterplot.matrix}{default}(x, labels=colnames(x), 
    diagonal=c("density", "boxplot", "histogram", "oned", "qqplot", "none"), 
    adjust=1, nclass, plot.points=TRUE, smooth=TRUE, span=0.5, reg.line=lm, 
    transform=FALSE, ellipse=FALSE, levels=c(.5, .9), robust=FALSE,
    groups=FALSE, by.groups=FALSE, col=palette(), 
    pch=1:n.groups, lwd=1, lwd.smooth=lwd,
    cex=par("cex"), cex.axis=par("cex.axis"), cex.labels=NULL, 
    cex.main=par("cex.main"),
    legend.plot=length(levels(groups)) > 1, ...)

spm(x, ...)
}

\arguments{
  \item{x}{a data matrix, numeric data frame, or formula.}
  \item{formula}{a one-side ``model'' formula, of the form 
    \code{ ~ x1 + x2 + ... + xk} or \code{ ~ x1 + x2 + ... + xk | z} where \code{z}
    evaluates to a factor or other variable to divide the data into groups.}
  \item{data}{for \code{scatterplot.matrix.formula},
    a data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{labels}{variable labels (for the diagonal of the plot).}
  \item{diagonal}{contents of the diagonal panels of the plot.}
  \item{adjust}{relative bandwidth for density estimate, passed to 
    \code{density} function.}
  \item{nclass}{number of bins for histogram, passed to \code{hist}
    function.}
  \item{plot.points}{if \code{TRUE} the points are plotted in each 
    off-diagonal panel.}
  \item{smooth}{if \code{TRUE} a lowess smooth is plotted in each
    off-diagonal panel.}
  \item{span}{span for lowess smoother.}
  \item{reg.line}{if not \code{FALSE} a line is plotted using the
    function given by this argument; e.g., using \code{rlm} in
    package \code{MASS} plots a robust-regression line.}
  \item{transform}{if \code{TRUE}, multivariate normalizing Box-Cox transformations
    are computed and plotted; if a vector of powers, one for each variable, these are applied
    as Box-Cox power transformations prior to plotting.}
  \item{ellipse}{if \code{TRUE} data-concentration ellipses are plotted in
    the off-diagonal panels.}
  \item{levels}{levels or levels at which concentration ellipses are plotted;
    the default is \code{c(.5, .9)}.}
  \item{robust}{if \code{TRUE} use the \code{cov.trob} function in the \code{MASS} package
    to calculate the center and covariance matrix for the data ellipse.}
  \item{groups}{a factor or other variable dividing the data into groups; groups are
    plotted with different colors and plotting characters.}
  \item{by.groups}{if \code{TRUE}, regression lines are fit by groups.}
  \item{pch}{plotting characters for points; default is the plotting characters in
    order (see \code{\link[graphics]{par}}).}
  \item{col}{colors for points and lines; the default is the
    in the current color palette, starting at the \emph{second} entry
    (see \code{\link[grDevices]{palette}} and \code{\link[graphics]{par}}).}
  \item{lwd}{width of linear-regression lines (default \code{1}).}
  \item{lwd.smooth}{width for smooth regression lines (default is the same as \code{lwd}).}
  \item{cex, cex.axis, cex.labels, cex.main}{set sizes of various graphical elements; 
    (see \code{\link[graphics]{par}}).}
  \item{legend.plot}{if \code{TRUE} then a legend for the groups is plotted
    in the bottom-right cell.}
  \item{...}{arguments to pass down.}
}

\value{
  \code{NULL}. This function is used for its side effect: producing
  a plot.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[graphics]{pairs}}, \code{\link{scatterplot}},
  \code{\link{data.ellipse}}, \code{\link{box.cox.powers}}, 
  \code{\link{box.cox}}, \code{\link[MASS]{cov.trob}}.}

\examples{
scatterplot.matrix(~income + education + prestige | type, data=Duncan)
scatterplot.matrix(~income + education + prestige, 
    transform=TRUE, data=Duncan)
}

\keyword{hplot}
