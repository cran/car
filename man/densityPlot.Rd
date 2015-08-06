\name{densityPlot}
\alias{densityPlot}
\alias{densityPlot.default}
\alias{densityPlot.formula}

\title{
Nonparametric Density Estimates
}
\description{
\code{densityPlot} contructs and graphs nonparametric density estimates, possibly conditioned on a factor.
It is a wrapper for the standard \R{} \code{\link{density}} function.
}

\usage{
densityPlot(x, ...)

\method{densityPlot}{default}(x, g, bw = "SJ", adjust=1,
    kernel = c("gaussian", "epanechnikov", "rectangular", 
               "triangular", "biweight", "cosine", "optcosine"),
    xlab = deparse(substitute(x)), ylab = "Density", col = palette(), 
    lty = seq_along(col), lwd = 2, grid=TRUE,
    legend.location = "topright", legend.title = deparse(substitute(g)), 
    show.bw = FALSE, rug = TRUE, ...)
    
\method{densityPlot}{formula}(formula, data = NULL, subset, na.action = NULL, xlab, ylab, ...)
}


\arguments{
  \item{x}{a numeric variable, the density of which is estimated.}
  \item{g}{an optional factor to divide the data.}
  \item{formula}{an \R{} model formula, of the form \code{~ variable} to estimate the unconditional
    density of \code{variable}, or \code{variable ~ factor} to estimate the density of \code{variable}
    within each level of \code{factor}.}
  \item{data}{an optional data frame containing the data.}
  \item{subset}{an optional vector defining a subset of the data.}
  \item{na.action}{a function to handle missing values; defaults to the value of the \R{} \code{na.action} option,
    initially set to \code{\link{na.omit}}.}
  \item{bw}{the bandwidth for the density estimate(s); either the quoted name of a rule to
    compute the bandwidth, or a numeric value; the default is \code{"SJ"}; if plotting by groups, \code{bw}
    may be a vector of values or rules, one for each group. See \code{\link{density}} and \code{\link{bw.SJ}} for details.}
  \item{adjust}{a multiplicative adjustment factor for the bandwidth; the default, \code{1}, indicates no adjustment;
    if plotting by groups, \code{adjust} may be a vector of adjustment factors, one for each group.}
  \item{kernel}{kernel function; the default is "gaussian" (see \code{\link{density}}).}
  \item{xlab}{label for the horizontal-axis; defaults to the name of the variable \code{x}.}
  \item{ylab}{label for the vertical axis; defaults to \code{"Density"}.}
  \item{col}{vector of colors for the density estimate(s); defaults to the color \code{\link{palette}}.}
  \item{lty}{vector of line types for the density estimate(s); defaults to the successive integers, starting at 1.}
  \item{lwd}{line width for the density estimate(s); defaults to 2.}
  \item{grid}{if \code{TRUE} (the default), grid lines are drawn on the plot.}
  \item{legend.location}{location for the legend when densities are plotted for several groups; defaults to
    \code{"upperright"}; see \code{\link{legend}}.}
  \item{legend.title}{label for the legend, which is drawn if densities are plotted by groups; the default is
    the name of the factor {g}.}
  \item{show.bw}{if \code{TRUE}, show the bandwidth(s) in the horizontal-axis label or (for multiple groups)
    the legend; the default is \code{FALSE}.}
  \item{rug}{if \code{TRUE} (the default), draw a rug plot (one-dimentional scatterplot) at the bottom of the density estimate.}
  \item{\dots}{arguments to be passed to \code{\link{plot}}.}
}

\value{
These functions return \code{NULL} invisibly and draw graphs.
}

\references{
W. N. Venables and B. D. Ripley (2002) \emph{Modern Applied Statistics with S}. New York: Springer.}

\author{
John Fox \email{jfox@mcmaster.ca}
}

\seealso{
\code{\link{density}}, \code{\link{bw.SJ}}, \code{link{plot.density}}
}

\examples{
densityPlot(~ income, show.bw=TRUE, data=Prestige)
densityPlot(income ~ type, data=Prestige)
}

\keyword{hplot}
