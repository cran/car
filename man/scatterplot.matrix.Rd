\name{scatterplot.matrix}
\alias{scatterplot.matrix}
\alias{spm}
\title{Scatterplot Matrices}
\description{
  Scatterplot matrices with univariate displays down the diagonal;
  \code{spm} is an abbreviation for \code{scatterplot.matrix}.
  This function just sets up a call to \code{pairs}.
}
\usage{

scatterplot.matrix(x, ...)

scatterplot.matrix.formula(formula, data=NULL, subset, ...)

scatterplot.matrix.default(data, labels=colnames(data), 
diagonal=c("density", "boxplot", "histogram", "qqplot", "none"), 
adjust=1, nclass, plot.points=TRUE, smooth=TRUE, span=0.5, reg.line=lm, 
transform=FALSE, groups=FALSE, by.groups=FALSE, col=palette(), 
pch=1:n.groups, lwd=1, legend.plot=length(levels(groups)) > 1)

spm(x, ...)
}

\arguments{
  \item{formula}{a one-side ``model'' formula, of the form 
    \code{ ~ x1 + x2 + ... + xk} or \code{ ~ x1 + x2 + ... + xk | z} where \code{z}
    evaluates to a factor or other variable to divide the data into groups.}
  \item{data}{data matrix or numeric data frame, or (for \code{scatterplot.matrix.formula})
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
    order (see \code{\link[base]{par}}).}
  \item{col}{colors for points and lines; the default is the
    in the current color palette, starting at the \emph{second} entry
    (see \code{\link[base]{palette}} and \code{\link[base]{par}}).}
  \item{lwd}{width for lines.}
  \item{legend.plot}{if \code{TRUE} then a legend for the groups is plotted
    in the bottom-right cell.}
}

\value{
  \code{NULL}. This function is used for its side effect: producing
  a plot.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[base]{pairs}}, \code{\link{scatterplot}},
  \code{\link{data.ellipse}}, \code{\link{box.cox.powers}}, 
  \code{\link{box.cox}}, \code{\link[MASS]{cov.trob}}.}

\examples{
data(Duncan)
scatterplot.matrix(~income + education + prestige | type, data=Duncan)
scatterplot.matrix(~income + education + prestige, 
    transform=T, data=Duncan)
}

\keyword{hplot}
