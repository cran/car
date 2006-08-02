\name{ceres.plots}
\alias{ceres.plots}
\alias{ceres.plot}
\alias{ceres.plot.lm}
\alias{ceres.plot.glm}

\title{Ceres Plots}
\description{
  These functions calculate Ceres plots for linear and generalized
  linear model.
}
\usage{
ceres.plots(model, variable, ask=missing(variable), one.page=!ask, 
  span=0.5, ...)

ceres.plot(model, ...)

\method{ceres.plot}{lm}(model, variable, line=TRUE, smooth=TRUE, span=0.5, iter,
  las=par('las'), col=palette()[2], pch=1, lwd=2, main="Ceres Plot", ...)

\method{ceres.plot}{glm}(model, ...)
}

\arguments{
  \item{model}{model object produced by \code{lm} or \code{glm}.}
  \item{variable}{variable (if it exists in the search path) or
    name of variable. This argument usually is omitted for
    \code{ceres.plots}.}
  \item{ask}{if \code{TRUE}, a menu is provided in the R Console for the
    user to select the variable(s) to plot, and to modify the span for the smoother
    used to draw a nonparametric-regression line on the plot.}
  \item{one.page}{if \code{TRUE} (and \code{ask=FALSE}), put all plots on one
    graph.}
  \item{span}{span for lowess smoother.}
  \item{iter}{number of robustness iterations for nonparametric-regression smooth;
    defaults to 3 for a linear model and to 0 for a non-Gaussian glm.}
  \item{line}{\code{TRUE} to plot least-squares line.}
  \item{smooth}{\code{TRUE} to plot nonparametric-regression (lowess) line.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{col}{color for points and lines; the default is the \emph{second} entry
    in the current color palette (see \code{\link[grDevices]{palette}}
    and \code{\link[graphics]{par}}).}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{main}{title for plot.}
  \item{\dots}{pass arguments down.}
}
\details{
  Ceres plots are a generalization of component+residual (partial
  residual) plots that are less prone to leakage of nonlinearity
  among the predictors.
  
  The function intended for direct use is \code{ceres.plots}. 
  By default, this function is used interactively
  through a text menu.
  
  The model cannot contain interactions, but can contain factors. 
  Factors may be present in the model, but Ceres plots cannot be drawn
  for them.
}
\value{
  \code{NULL}. These functions are used for their side effect: producing
  plots.
}
\references{ 
  Cook, R. D. and Weisberg, S. (1999)
  \emph{Applied Regression, Including Computing and Graphics.} Wiley.

  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{cr.plots}}, \code{\link{av.plots}}}

\examples{
  \dontrun{
ceres.plots(lm(prestige~income+education+type, data=Prestige))
  }
}

\keyword{hplot}
\keyword{regression}
