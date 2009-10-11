\name{cr.plots}
\alias{cr.plots}
\alias{crp}
\alias{cr.plot}
\alias{cr.plot.lm}
\alias{cr.plot.glm}

\title{Component+Residual (Partial Residual) Plots}
\description{
  These functions construct component+residual plots (also called
  partial-residual plots) for linear and generalized linear models.
}
\usage{
cr.plots(model, variable, ask=missing(variable), one.page=!ask, 
  span=0.5, ...)

crp(...)

cr.plot(model, ...)

\method{cr.plot}{lm}(model, variable, order=1, line=TRUE, smooth=TRUE, iter, 
  span=0.5, las=par('las'), col=palette()[2], pch=1, lwd=2,
  main="Component+Residual Plot", ...)

\method{cr.plot}{glm}(model, ...)
}

\arguments{
  \item{model}{model object produced by \code{lm} or \code{glm}.}
  \item{variable}{variable (if it exists in the search path) or
    name of variable. This argument usually is omitted for \code{crp}
    or \code{cr.plots}.}
  \item{ask}{if \code{TRUE}, a menu is provided in the R Console for the
    user to select the variable(s) to plot, and to modify the span for the smoother
    used to draw a nonparametric-regression line on the plot.}
  \item{one.page}{if \code{TRUE} (and \code{ask=FALSE}), put all plots on one
    graph.}
  \item{order}{order of polynomial regression performed for predictor to be plotted.}
  \item{line}{\code{TRUE} to plot least-squares line.}
  \item{smooth}{\code{TRUE} to plot nonparametric-regression (lowess) line.}
  \item{iter}{number of robustness iterations for nonparametric-regression smooth;
    defaults to 3 for a linear model and to 0 for a non-Gaussian glm.}
  \item{span}{span for lowess smoother.}
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
  The function intended for direct use is \code{cr.plots} (for which \code{crp}
  is an abbreviation). By default, these functions are used interactively
  through a text menu.
  
  The model cannot contain interactions, but can contain factors. 
  Parallel boxplots of the partial residuals are drawn for the levels
  of a factor.
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

\seealso{\code{\link{ceres.plots}}, \code{\link{av.plots}}}

\examples{
cr.plots(lm(prestige~income+education, data=Prestige), 
  variable="income")
  \dontrun{
cr.plots(glm(partic != "not.work" ~ hincome + children, 
  data=Womenlf, family=binomial))
  }
}

\keyword{hplot}
\keyword{regression}
