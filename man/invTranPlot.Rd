\name{invTranPlot}
\alias{invTranPlot}
\alias{invTranPlot.default}
\alias{invTranPlot.formula}
\alias{invTranEstimate}

\title{ Choose a Predictor Transformation Visually or Numerically }
\description{
\code{invTranPlot}
draws a two-dimensional scatterplot of \eqn{Y}{Y} versus 
\eqn{X}{X}, along with the OLS
fit from the regression of \eqn{Y}{Y} on 
\eqn{(X^{\lambda}-1)/\lambda}{(Y^(lam)-1)/lam}.  \code{invTranEstimate}
finds the nonlinear least squares estimate of \eqn{\lambda}{lambda} and its
standard error.
}
\usage{            
invTranPlot(x, ...)

\S3method{invTranPlot}{formula}(x, data, subset, na.action, ...)

\S3method{invTranPlot}{default}(x, y, lambda=c(-1, 0, 1), robust=FALSE,
        lty.lines=rep(c("solid", "dashed", "dotdash", "longdash", "twodash"), 
        length=1 + length(lambda)), lwd.lines=2, 
        col=palette()[1], col.lines=palette(),
        xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
        family="bcPower", optimal=TRUE, key="auto",
        id.method = "x",
        labels, 
        id.n = if(id.method[1]=="identify") Inf else 0,
        id.cex=1, id.col=palette()[1], grid=TRUE, ...)

invTranEstimate(x, y, family="bcPower", confidence=0.95, robust=FALSE)
}

\arguments{
  \item{x}{The predictor variable, or a formula with a single response and
  a single predictor }
  \item{y}{The response variable }
  \item{data}{An optional data frame to get the data for the formula}
  \item{subset}{Optional, as in \code{lm}, select a subset of the cases}
  \item{na.action}{Optional, as in \code{lm}, the action for missing data}
  \item{lambda}{The powers used in the plot.  The optimal power than minimizes
  the residual sum of squares is always added unless optimal is \code{FALSE}. }
  \item{robust}{If TRUE, then the estimated transformation is computed using
Huber M-estimation with the MAD used to estimate scale and k=1.345.  The
default is FALSE.} 
  \item{family}{The transformation family to use, \code{"bcPower"}, 
  \code{"yjPower"}, or a user-defined family.}
  \item{confidence}{returns a profile likelihood confidence interval for the optimal 
	transformation with this confidence level.  If \code{FALSE}, or if \code{robust=TRUE},
no interval is returned.}
  \item{optimal}{Include the optimal value of lambda?}
  \item{lty.lines}{line types corresponding to the powers}
  \item{lwd.lines}{the width of the plotted lines, defaults to 2 times the standard}
  \item{col}{color(s) of the points in the plot.  If you wish to distinguish points
  according to the levels of a factor, we recommend using symbols, specified with
  the \code{pch} argument, rather than colors.}
  \item{col.lines}{color of the fitted lines corresponding to the powers.  The
  default is to use the colors returned by  \code{\link{palette}}}
  \item{key}{The default is \code{"auto"}, in which case a legend is added to
the plot, either above the top marign or in the bottom right or top right corner.
Set to NULL to suppress the legend.}
  \item{xlab}{Label for the horizontal axis.}

  \item{ylab}{Label for the vertical axis.}
  \item{id.method,labels,id.n,id.cex,id.col}{Arguments for the labelling of
    points.  The default is \code{id.n=0} for labeling no points.  See
    \code{\link{showLabels}} for details of these arguments.
}
  \item{...}{Additional arguments passed to the plot method, such as \code{pch}.}
\item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
}
\value{
\code{invTranPlot}
  plots a graph and returns a data frame with \eqn{\lambda}{lam} in the 
  first column, and the residual sum of squares from the regression
  for that \eqn{\lambda}{lam} in the second column.

  \code{invTranEstimate} returns a list with elements \code{lambda} for the
  estimate, \code{se} for its standard error, and \code{RSS}, the minimum
  value of the residual sum of squares.  
}
\seealso{ \code{\link{inverseResponsePlot}},\code{\link{optimize}}}

\references{
Fox, J. and Weisberg, S. (2011) 
\emph{An R Companion to Applied Regression}, Second Edition, Sage.

Pendergast, L, and Sheather, S. (in press).  On sensitivity of response plot
estimation of a robust estimation approach.  \emph{Scandinavian Journal of
Statistics}.

Weisberg, S. (2014)  \emph{Applied Linear Regression}, Fourth Edition, Wiley
  Wiley. 
}

\author{Sanford Weisberg, \email{sandy@umn.edu} }

\examples{
with(UN, invTranPlot(gdp, infant.mortality))
with(UN, invTranEstimate(gdp, infant.mortality))
}    
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{regression}
