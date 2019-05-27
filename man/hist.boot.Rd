\name{hist.boot}
\alias{hist.boot}
\alias{summary.boot}
\alias{confint.boot}
\alias{vcov.boot}
\alias{Confint.boot}

\title{
Methods Functions to Support \code{boot} Objects
}
\description{
The \code{Boot} function in the \pkg{car} package uses the \code{\link{boot}} function from the
\pkg{boot} package to do a straightforward case
or residual bootstrap for many regression objects.  These are method functions for standard generics to summarize the results of the bootstrap.  Other tools for this purpose are available in the \code{boot} package.
}

\usage{
\method{hist}{boot}(x, parm, layout = NULL, ask, main = "", freq = FALSE,
    estPoint = TRUE, point.col = carPalette()[1], point.lty = 2, point.lwd = 2,
    estDensity = !freq, den.col = carPalette()[2], den.lty = 1, den.lwd = 2,
    estNormal = !freq, nor.col = carPalette()[3], nor.lty = 2, nor.lwd = 2,
    ci = c("bca", "none", "perc", "norm"), level = 0.95, 
    legend = c("top", "none", "separate"), box = TRUE, ...)

\method{summary}{boot}(object, parm, high.moments = FALSE, extremes = FALSE, ...)

\method{confint}{boot}(object, parm, level = 0.95, type = c("bca", "norm",
    "basic", "perc"), ...)

\method{Confint}{boot}(object, parm, level = 0.95, type = c("bca", "norm",
    "basic", "perc"), ...)

\method{vcov}{boot}(object, ...)
}

\arguments{
  \item{x, object}{
An object created by a call to \code{boot} in the \code{boot} package, or to \code{Boot} in the \pkg{car} package of class \code{"boot"}.
}
  \item{parm}{
A vector of numbers or coefficient names giving the coefficients for which a
histogram or confidence interval is desired.  If numbers are used, 1 corresponds to the intercept,
if any.  The default is all coefficients.
}
  \item{layout}{
If set to a value like \code{c(1, 1)} or \code{c(4, 3)}, the layout
of the graph will have this many rows and columns.  If not set, the program
will select an appropriate layout.  If the number of graphs exceed nine, you
must select the layout yourself, or you will get a maximum of nine per page.
If \code{layout=NA}, the function does not set the layout and the user can
use the \code{par} function to control the layout, for example to have
plots from two models in the same graphics window.
}
  \item{ask}{
If \code{TRUE}, ask the user before drawing the next plot; if \code{FALSE}, don't
ask.
}
  \item{main}{
Main title for the graphs.  The default is \code{main=""} for no title.
}
  \item{freq}{
The default for the generic \code{hist} function is \code{freq=TRUE} to give a frequency histogram.  The default for \code{hist.boot} is \code{freq=FALSE} to give a density histogram.  A density estimate and/or a fitted normal density can be added to the graph if \code{freq=FALSE} but not if \code{freq=TRUE}.
}
  \item{estPoint, point.col, point.lty, point.lwd}{
If \code{estPoint=TRUE}, the default, a vertical line is drawn on the histgram
at the value of the point estimate computed from the complete data. The
remaining three optional arguments set the color, line type and line width
of the line that is drawn.
}

  \item{estDensity, den.col, den.lty, den.lwd}{
If \code{estDensity=TRUE} and\code{freq=FALSE}, the default, a kernel density estimate is drawn
on the plot with a call to the \code{density} function with no additional
arguments. The
remaining three optional arguments set the color, line type and line width
of the lines that are drawn.
}
  \item{estNormal, nor.col, nor.lty, nor.lwd}{
If \code{estNormal=TRUE} and\code{freq=FALSE}, the default, a normal density
with mean and sd computed from the data is drawn on the plot. The
remaining three optional arguments set the color, line type and line width
of the lines that are drawn.
}
  \item{ci}{
A confidence interval based on the bootstrap will be added to the histogram
using the BCa method if \code{ci="bca"} the percentile method if
\code{ci="perc"}, or the normal method if \code{ci="norm"}.  No interval is drawn if
\code{ci="none"}.  The default is \code{"bca"}.  The interval is indicated
by a thick horizontal line at \code{y=0}.  For some bootstraps the BCa method is unavailable, in which case a warning is issued and \code{ci="perc"} is substituted.  If you wish to see all the options at once, see \code{\link[boot]{boot.ci}}.  The normal method is computed as the (estimate from the original data) minus the bootstrap bias plus or minus the standard deviation of the bootstrap replicates times the appropriate quantile of the standard normal distribution.
}
  \item{legend}{
A legend can be added to the (array of) histograms.  The value \dQuote{"top"}
puts at the top-left of the plots.  The value \dQuote{"separate"} puts the legend in its own graph following all the histograms.  The value \dQuote{"none"}
suppresses the legend.
}
  \item{box}{
Add a box around each histogram.
}
  \item{\dots}{
Additional arguments passed to \code{hist}; for other methods this is included for compatibility with the generic method.  For example, the argument
\code{border=par()$bg} in \code{hist} will draw the histogram transparently, leaving only
the density estimates.  With the \code{vcov} function, the additional arguments are passed to \code{\link{cov}}.  See the Value section, below.
}
  \item{high.moments}{
Should the skewness and kurtosis be included in the summary?  Default is FALSE.
}
  \item{extremes}{
Should the minimum, maximum and range be included in the summary?  Default is FALSE.
}
  \item{level}{
Confidence level, a number between 0 and 1.  In \code{confint}, \code{level} can be
a vector; for example \code{level=c(.50, .90, .95)} will return the following estimated quantiles:  \code{c(.025, .05, .25, .75, .95, .975)}.
}
  \item{type}{
Selects the  confidence interval type.  The types
implemented are the \code{"percentile"} method, which uses the function
\code{quantile} to return the appropriate quantiles for the confidence
limit specified, the default \code{bca} which uses the bias-corrected and accelerated
method presented by Efron and Tibshirani (1993, Chapter 14).  For the
other types, see the documentation for \code{\link[boot]{boot}}.
}
}
\value{
\code{hist} is used for the side-effect of drawing an array of historgams of
each column of the first argument.   \code{summary} returns a matrix of
summary statistics for each of the columns in the bootstrap object.  The
\code{confint} method returns confidence intervals.  \code{Confint} appends the estimates based on the original fitted model to the left of the confidence intervals.

The \code{vcov} returns the sample covariance of the bootstrap sample estimates.  If any of the bootstrap replications returned an \code{NA} value, then \code{vcov} will return a matrix of \code{NA}s.  You can compute the covariance of the complete cases by adding the argument  \code{use="complete.obs"} to the call to \code{vcov}; see \code{\link{cov}}.
}
\references{

Efron, B. and Tibsharini, R. (1993)
\emph{An Introduction to the Bootstrap}. New
York: Chapman and Hall.

Fox, J. and Weisberg, S. (2019)
\emph{An R Companion to Applied Regression}, Third Edition.  Thousand Oaks:  Sage.

Fox, J. and Weisberg, S. (2018) \emph{Bootstrapping Regression Models in R},
\url{https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Bootstrapping.pdf}.

Weisberg, S. (2013) \emph{Applied
Linear Regression}, Fourth Edition, Wiley
}
\author{Sanford Weisberg, \email{sandy@umn.edu}}

\seealso{See Also \code{\link{Boot}}, \code{\link{hist}},
\code{\link{density}}, Fox and Weisberg (2017), cited above }

\examples{
m1 <- lm(Fertility ~ ., swiss)
betahat.boot <- Boot(m1, R=99) # 99 bootstrap samples--too small to be useful
summary(betahat.boot)  # default summary
confint(betahat.boot)
hist(betahat.boot)
}
\keyword{regression}

