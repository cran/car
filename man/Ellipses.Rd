\name{Ellipses}
\alias{ellipse}
\alias{dataEllipse}
\alias{confidenceEllipse}
\alias{confidenceEllipse.lm}
\alias{confidenceEllipse.glm}

\title{Ellipses, Data Ellipses, and Confidence Ellipses}
\description{
  These functions draw ellipses, including data ellipses, and
  confidence ellipses for linear and generalized linear models.
}
\usage{
ellipse(center, shape, radius, log="", center.pch=19, center.cex=1.5, 
  segments=51, draw=TRUE, add=draw, xlab="", ylab="", 
   col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, grid=TRUE, ...)

dataEllipse(x, y, weights, log="", levels=c(0.5, 0.95), center.pch=19, center.cex=1.5, 
  draw=TRUE, plot.points=draw, add=!plot.points, segments=51, robust=FALSE,
  xlab=deparse(substitute(x)), 
  ylab=deparse(substitute(y)), 
  col=palette()[1:2], lwd=2, fill=FALSE, fill.alpha=0.3, grid=TRUE, ...)

confidenceEllipse(model, ...)

\method{confidenceEllipse}{lm}(model, which.coef, levels=0.95, Scheffe=FALSE, dfn,
  center.pch=19, center.cex=1.5, segments=51, xlab, ylab, 
  col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, draw=TRUE, add=!draw, ...)

\method{confidenceEllipse}{glm}(model, which.coef, levels=0.95, Scheffe=FALSE,  dfn,
  center.pch=19, center.cex=1.5, segments=51, xlab, ylab, 
  col=palette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, draw=TRUE, add=!draw, ...)
}

\arguments{
  \item{center}{2-element vector with coordinates of center of ellipse.}
  \item{shape}{\eqn{2\times 2}{2 * 2} shape (or covariance) matrix.}
  \item{radius}{radius of circle generating the ellipse.}
  \item{log}{when an ellipse is to be added to an existing plot, indicates
    whether computations were on logged values and to be plotted on logged
    axes; \code{"x"} if the x-axis is logged, \code{"y"} if the y-axis is
    logged, and \code{"xy"} or \code{"yx"} if both axes are logged. The
    default is \code{""}, indicating that neither axis is logged.}
  \item{center.pch}{character for plotting ellipse center.}
  \item{center.cex}{relative size of character for plotting ellipse center.}
  \item{segments}{number of line-segments used to draw ellipse.}
  \item{draw}{if \code{TRUE} produce graphical output; if \code{FALSE}, only invisibly return coordinates of ellipse(s).}
  \item{add}{if \code{TRUE} add ellipse(s) to current plot.}
  \item{xlab}{label for horizontal axis.}
  \item{ylab}{label for vertical axis.}
  \item{x}{a numeric vector, or (if \code{y} is missing) a 2-column numeric matrix.}
  \item{y}{a numeric vector, of the same length as \code{x}.}
  \item{weights}{a numeric vector of weights, of the same length as \code{x} and \code{y} to be used  by \code{\link[stats]{cov.wt}}
    or \code{\link[MASS]{cov.trob}} in computing a weighted covariance matrix; if absent, weights of \code{1} are used.}
  \item{plot.points}{if \code{FALSE} data ellipses are drawn,
    but points are not plotted.}
  \item{levels}{draw elliptical contours at these (normal) probability or confidence levels.}
  \item{robust}{if \code{TRUE} use the \code{cov.trob} function in the \pkg{MASS} package
    to calculate the center and covariance matrix for the data ellipse.}
  \item{model}{a model object produced by \code{lm} or \code{glm}.}
  \item{which.coef}{2-element vector giving indices of coefficients to plot; 
    if missing, the first two coefficients (disregarding the regression constant)
    will be selected.}
  \item{Scheffe}{if \code{TRUE} scale the ellipse so that its projections onto the
    axes give Scheffe confidence intervals for the coefficients.}
  \item{dfn}{``numerator'' degrees of freedom (or just degrees of freedom for a GLM) for
  	drawing the confidence ellipse. Defaults to the number of coefficients in the model (disregarding the constant) if 
  	\code{Scheffe} is \code{TRUE}, or to \code{2} otherwise; selecting \code{dfn = 1} will 
  	draw the ``confidence-interval generating'' ellipse, with projections on the axes
  	corresponding to individual confidence intervals with the stated level of coverage.}
  \item{col}{color for lines and ellipse center; the default is the \emph{second} entry
    in the current color palette (see \code{\link[grDevices]{palette}}
    and \code{\link[graphics]{par}}). For \code{dataEllipse}, two colors can be given, in
    which case the first is for plotted points and the second for lines and the ellipse center.}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{fill}{fill the ellipse with translucent color \code{col} (default, \code{FALSE})?}
  \item{fill.alpha}{transparency of fill (default = \code{0.3}).}
  \item{\dots}{other plotting parameters to be passed to \code{plot} and
    \code{line}.}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
}
\details{
  The ellipse is computed by suitably transforming a unit circle.
  
  \code{dataEllipse} superimposes the normal-probability contours over a scatterplot
  of the data.
}
\value{
   These functions are mainly used for their side effect of producing plots.  For
   greater flexibility (e.g., adding plot annotations), however, \code{ellipse} returns invisibly the (x, y) coordinates of the calculated ellipse.
   \code{dataEllipse} and \code{confidenceEllipse} return invisibly the  coordinates of one or more ellipses, in the latter instance a list named by
   \code{levels}.
}

\references{ 
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Monette, G. (1990)
  Geometry of multiple regression and 3D graphics.
  In Fox, J. and Long, J. S. (Eds.) 
  \emph{Modern Methods of Data Analysis.} Sage.
 }
 
\author{Georges Monette, John Fox \email{jfox@mcmaster.ca}, and Michael Friendly.}

\seealso{\code{\link[MASS]{cov.trob}}, \code{\link[stats]{cov.wt}}.}
  
\examples{
dataEllipse(Prestige$income, Prestige$education, levels=0.1*1:9, lty=2, 
	fill=TRUE, fill.alpha=0.1)
confidenceEllipse(lm(prestige~income+education, data=Prestige), Scheffe=TRUE)

wts <- rep(1, nrow(Duncan))
wts[c(6, 16)] <- 0 # delete Minister, Conductor
with(Duncan, {
	dataEllipse(income, prestige, levels=0.68)
	dataEllipse(income, prestige, levels=0.68, robust=TRUE, plot.points=FALSE, col="green3")
	dataEllipse(income, prestige, weights=wts, levels=0.68, plot.points=FALSE, col="brown")
	dataEllipse(income, prestige, weights=wts, robust=TRUE, levels=0.68, 
		plot.points=FALSE, col="blue")
	})
}

\keyword{hplot}
\keyword{aplot}
