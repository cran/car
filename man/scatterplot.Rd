\name{scatterplot}
\alias{scatterplot}
\alias{scatterplot.formula}
\alias{scatterplot.default}
\alias{sp}

\title{Scatterplots with Boxplots}
\description{
  Makes fancy scatterplots, with boxplots in the margins, a lowess
  smooth, and a regression line; \code{sp} is an abbreviation for
  \code{scatterplot}.
}
\usage{
scatterplot(x, ...)

scatterplot.formula(formula, data, xlab, ylab, subset, labels=FALSE, ...)

scatterplot.default(x, y, smooth=TRUE, span=0.5, reg.line=lm, 
  boxplots="xy", xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), 
  las=par('las'), lwd=1, labels=FALSE, log="", groups=FALSE, by.groups=!(groups == FALSE),
  ellipse=FALSE, levels=c(.5, .9), robust=FALSE, 
  col=palette(), pch=1:n.groups, legend.plot=length(levels(groups)) > 1, 
  reset.par=TRUE, ...)
  
sp{...}
}

\arguments{
  \item{formula}{``model'' formula, of the form \code{y ~ x} or 
    (to plot by groups) \code{y ~ x | z}, where \code{z} evaluates to a factor 
    or other variable dividing the data into groups.}
  \item{data}{data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{x}{vector of horizontal coordinates.}
  \item{y}{vector of verical coordinates.}
  \item{smooth}{if \code{TRUE} a lowess nonparametric regression line is
    drawn on the plot.}
  \item{span}{span for the lowess smooth.}
  \item{reg.line}{function to draw a regression line on the plot
    or \code{FALSE} not to plot a regression line.}
  \item{boxplots}{if \code{"x"} a boxplot for \code{x} is drawn above the plot;
    if \code{"y"} a boxplot for \code{y} is drawn to the right of the plot;
    if \code{"xy"} both boxplots are drawn.}
  \item{xlab}{label for horizontal axis.}
  \item{ylab}{label for vertical axis.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[base]{par}}).}
  \item{lwd}{width of plotted lines.}
  \item{labels}{if not \code{FALSE} a vector of point labels, to be used interactively
    to identify points on the plot.}
  \item{log}{same as the \code{log} argument to \code{plot}, to produce log axes.}
  \item{groups}{a factor or other variable dividing the data into groups; groups are
    plotted with different colors and plotting characters.}
  \item{by.groups}{if \code{TRUE}, regression lines are fit by groups.}
  \item{ellipse}{if \code{TRUE} data-concentration ellipses are plotted.}
  \item{levels}{level or levels at which concentration ellipses are plotted;
    the default is \code{c(.5, .9)}.}
  \item{robust}{if \code{TRUE} use the \code{cov.trob} function in the \code{MASS} package
    to calculate the center and covariance matrix for the data ellipse.}
  \item{col}{colors for points and lines; the default is the
    current color palette, starting at the \emph{second} entry
    (see \code{\link[base]{palette}} and \code{\link[base]{par}}).}
  \item{pch}{plotting characters for points; default is the plotting characters in
    order (see \code{\link[base]{par}}).}
  \item{legend.plot}{if \code{TRUE} then a legend for the groups is plotted, and positioned
    interactively with a mouse click; if points are to be labelled, then the legend is 
    positioned first.}
  \item{reset.par}{if \code{TRUE} then plotting parameters are reset to their previous values
    when \code{scatterplot} exits; if \code{FALSE} then the \code{mar} and \code{mfcol} parameters are
    altered for the current plotting device. Set to \code{FALSE} if you want to add graphical elements
    (such as lines) to the plot.}
  \item{\dots}{other arguments passed to \code{plot}.}
}

\value{
  \code{NULL}. This function is used for its side effect: producing
  a plot.
}

\author{John Fox \email{jfox@mcmaster.ca}}


\seealso{\code{\link[base]{plot}}, \code{\link[base]{boxplot}}, 
  \code{\link{scatterplot.matrix}}, \code{\link{data.ellipse}},
  \code{\link[base]{par}}, \code{\link[MASS]{cov.trob}}.}

\examples{
data(Prestige)
scatterplot(prestige~income|type, data=Prestige, span=1, legend=FALSE)
  \dontrun{
data(UN)
scatterplot(infant.mortality~gdp, labels=row.names(UN), data=UN)

scatterplot(infant.mortality~gdp, log="xy", labels=row.names(UN), data=UN)
  }
}

\keyword{hplot}
