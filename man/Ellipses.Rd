\name{Ellipses}
\alias{ellipse}
\alias{data.ellipse}
\alias{confidence.ellipse}
\alias{confidence.ellipse.lm}
\alias{confidence.ellipse.glm}

\title{Ellipses, Data Ellipses, and Confidence Ellipses}
\description{
  These functions draw ellipses, including data ellipses, and
  confidence ellipses for linear and generalized linear models.
}
\usage{
ellipse(center, shape, radius, center.pch=19, center.cex=1.5, 
  segments=51, add=TRUE, xlab="", ylab="", 
   las=par('las'), col=palette()[2], lwd=2, lty=1, ...)

data.ellipse(x, y, levels=c(0.5, 0.9), center.pch=19, center.cex=1.5, 
  plot.points=TRUE, add=!plot.points, segments=51, robust=FALSE,
  xlab=deparse(substitute(x)), 
  ylab=deparse(substitute(y)), 
  las=par('las'), col=palette()[2], pch=1, lwd=2, lty=1, ...)

confidence.ellipse(model, ...)

\method{confidence.ellipse}{lm}(model, which.coef, levels=0.95, Scheffe=FALSE, 
  center.pch=19, center.cex=1.5, segments=51, xlab, ylab, 
  las=par('las'), col=palette()[2], lwd=2, lty=1, ...)

\method{confidence.ellipse}{glm}(model, which.coef, levels=0.95, Scheffe=FALSE, 
  center.pch=19, center.cex=1.5, segments=51, xlab, ylab, 
  las=par('las'), col=palette()[2], lwd=2, lty=1, ...)
}

\arguments{
  \item{center}{2-element vector with coordinates of center of ellipse.}
  \item{shape}{\eqn{2\times 2}{2 * 2} shape (or covariance) matrix.}
  \item{radius}{radius of circle generating the ellipse.}
  \item{center.pch}{character for plotting ellipse center.}
  \item{center.cex}{relative size of character for plotting ellipse center.}
  \item{segments}{number of line-segments used to draw ellipse.}
  \item{add}{if \code{TRUE} add ellipse to current plot.}
  \item{xlab}{label for horizontal axis.}
  \item{ylab}{label for vertical axis.}
  \item{x}{a numeric vector, or (if \code{y} is missing) a 2-column numeric matrix.}
  \item{y}{a numeric vector, of the same length as \code{x}.}
  \item{plot.points}{if \code{FALSE} data ellipses are added to the current scatterplot,
    but points are not plotted.}
  \item{levels}{draw elliptical contours at these (normal) probability or confidence levels.}
  \item{robust}{if \code{TRUE} use the \code{cov.trob} function in the \code{MASS} package
    to calculate the center and covariance matrix for the data ellipse.}
  \item{model}{a model object produced by \code{lm} or \code{glm}.}
  \item{which.coef}{2-element vector giving indices of coefficients to plot; 
    if missing, the first two coefficients (disregarding the regression constant)
    will be selected.}
  \item{Scheffe}{if \code{TRUE} scale the ellipse so that its projections onto the
    axes give Scheffe confidence intervals for the coefficients.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{col}{color for points and lines; the default is the \emph{second} entry
    in the current color palette (see \code{\link[grDevices]{palette}}
    and \code{\link[graphics]{par}}).}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{lty}{line type; default is \code{1}, a solid line (see \code{\link[graphics]{par}}).}
  \item{\dots}{other plotting parameters to be passed to \code{plot} and
    \code{line}.}
}
\details{
  The ellipse is computed by suitably transforming a unit circle.
  
  \code{data.ellipse} superimposes the normal-probability contours over a scatterplot
  of the data.
}
\value{
  \code{NULL}. These functions are used for their side effect: producing
  plots.
}
\references{ 
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
  
  Monette, G. (1990)
  Geometry of multiple regression and 3D graphics.
  In Fox, J. and Long, J. S. (Eds.) 
  \emph{Modern Methods of Data Analysis.} Sage.
 }
 
\author{Georges Monette \email{Georges.Monette@mathstat.YorkU.CA} 
  and John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[MASS]{cov.trob}}.}
  
\examples{
data.ellipse(Prestige$income, Prestige$education, levels=0.1*1:9, lty=2)
confidence.ellipse(lm(prestige~income+education, data=Prestige), Scheffe=TRUE)
}

\keyword{hplot}
\keyword{aplot}
