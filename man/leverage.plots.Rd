\name{leverage.plots}
\alias{leverage.plots}
\alias{leverage.plot}
\alias{leverage.plot.lm}
\alias{leverage.plot.glm}

\title{Regression Leverage Plots}
\description{
  These functions display a generalization, due to Sall (1990), 
  of added-variable plots to multiple-df terms in a linear model. 
  When a term has just 1 df, the leverage plot is a rescaled version
  of the usual added-variable (partial-regression) plot.
}
\usage{
leverage.plots(model, term.name, ask=missing(term.name), ...)

leverage.plot(model, ...)

\method{leverage.plot}{lm}(model, term.name, 
  labels=names(residuals(model)[!is.na(residuals(model))]),  
  identify.points=TRUE, las=par('las'), col=palette()[2], pch=1, lwd=2, 
  main="Leverage Plot", ...)

\method{leverage.plot}{glm}(model, ...)
}

\arguments{
  \item{model}{model object produced by \code{lm}}.
  \item{term.name}{name of term in the model to be plotted; this argument is usually
    omitted for \code{leverage.plots}.}
  \item{ask}{if \code{TRUE}, a menu is provided in the R Console for the
    user to select the term(s) to plot.}
  \item{labels}{observation names.}
  \item{identify.points}{if \code{TRUE}, then identify points interactively.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{col}{color for points and lines; the default is the \emph{second} entry
    in the current color palette (see \code{\link[grDevices]{palette}}
    and \code{\link[graphics]{par}}).}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{main}{title for plot.}
  \item{\dots}{arguments passed down to method functions.}
}
\details{
  The function intended for direct use is \code{leverage.plots}. By default, this
  function is used interactively through a text menu.
  
  The model can contain factors and interactions. A leverage plot can be
  drawn for each term in the model, including the constant.
  
  \code{leverage.plot.glm} is a dummy function, which generates an error message.
}

\value{
  \code{NULL}. These functions are used for their side effect: producing
  plots.
}

\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
  
  Sall, J. (1990)
  Leverage plots for general linear hypotheses.
  \emph{American Statistician} \bold{44}, 308--315.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{av.plots}}}

\examples{
  \dontrun{
leverage.plots(lm(prestige~(income+education)*type, data=Duncan))
  }
}

\keyword{hplot}
\keyword{regression}
