\name{invResPlot}
\alias{invResPlot}
\alias{inverseResponsePlot}
\alias{inverseResponsePlot.lm}

\title{Inverse Response Plots to Transform the Response}
\description{
For a \code{lm} model, draws an inverse.response plot with the response \eqn{Y}{Y} on the
vertical axis and the fitted values \eqn{\hat{Y}}{Yhat}
on the horizontal axis.  Uses \code{nls} to
estimate  \eqn{\lambda}{lambda} in the function 
\eqn{\hat{Y}=b_0+b_1Y^{\lambda}}{Yhat = b0 + b1(Y)^(lambda)}.
Adds the fitted curve to the plot. 
\code{invResPlot} is an alias for \code{inverseResponsePlot}.
}

\usage{
inverseResponsePlot(model, lambda=c(-1,0,1), robust=FALSE, xlab=NULL, ...)

\S3method{inverseResponsePlot}{lm}(model, lambda=c(-1,0,1), robust=FALSE, 
   xlab=NULL, labels=names(residuals(model)), ...)

invResPlot(model, ...)
}

\arguments{
  \item{model}{A \code{lm} regression object}
  \item{lambda}{A vector of values for lambda.  A plot will be produced 
  with curves corresponding to these lambdas and to the least squares 
  estimate of lambda}
  \item{xlab}{The horizontal axis label.  If \code{NULL}, it is constructed by the function.}
  \item{labels}{Case labels if labeling is turned on; see
\code{\link{invTranPlot}} and \code{\link{showLabels}} for arguments.}
  \item{robust}{If TRUE, then estimation uses Huber M-estimates with the median
absolute deviation to estimate scale and k= 1.345.  The default is FALSE.} 
  \item{\dots}{Other arguments passed to \code{invTranPlot} and then to
\code{plot}.}
}

\value{
As a side effect, a plot is produced with the response on the horizontal
axis and fitted values on the vertical axis.  
Several lines are added to be plot as
the ols estimates of the regression of \eqn{\hat{Y}}{Yhat} on
\eqn{Y^{\lambda}}{Y^(lambda)},
interpreting \eqn{\lambda}{lambda} = 0 to be natural logarithms.

Numeric output is a list with elements
  \item{lambda}{Estimate of transformation parameter for the response}
  \item{RSS}{The residual sum of squares at the minimum if robust=FALSE.
If robust = TRUE, the value of Huber objective function is returned.}
}
\seealso{\code{\link{invTranPlot}}, \code{\link{powerTransform}}, 
\code{\link{showLabels}}}

\references{
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Pendergast, L, and Sheather, S. (in press).  On sensitivity of response plot
estimation of a robust estimation approach.  \emph{Scandinavian Journal of
Statistics}.

  Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley, Chapter 7.
}

\author{Sanford Weisberg, \code{sandy@umn.edu}}

\examples{
m2 <- lm(rate ~ log(len) + log(ADT) + slim + shld + log(sigs1), Highway1)
invResPlot(m2)
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ regression}
