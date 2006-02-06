\name{av.plots}
\alias{av.plots}
\alias{avp}
\alias{av.plot}
\alias{av.plot.lm}
\alias{av.plot.glm}

\title{Added-Variable Plots}
\description{
  These functions construct added-variable (also called 
  partial-regression) plots for  linear and generalized linear
  models.
}

\usage{
av.plots(model, variable, ask=missing(variable), one.page=!ask, ...)

avp(...)

av.plot(model, ...)

\method{av.plot}{lm}(model, variable, 
  labels=names(residuals(model)[!is.na(residuals(model))]), 
  identify.points=TRUE, las=par("las"), col=palette()[2], pch=1, lwd=2,
  main="Added-Variable Plot", ...)

\method{av.plot}{glm}(model, variable, 
    labels=names(residuals(model)[!is.na(residuals(model))]), 
    identify.points=TRUE, las=par("las"), col=palette()[2], pch=1, lwd=2, 
    main="Added-Variable Plot", type=c("Wang", "Weisberg"), ...)
}

\arguments{
  \item{model}{model object produced by \code{lm} or \code{glm}.}
  \item{variable}{variable (if it exists in the search path) or
    name of variable. This argument usually is omitted for \code{avp}
    or \code{av.plots}.}
  \item{ask}{if \code{TRUE}, a menu is provided in the R Console for the
    user to select the term(s) to plot.}
  \item{one.page}{if \code{TRUE} (and \code{ask=FALSE}), put all plots on one
    graph.}
  \item{labels}{observation names.}
  \item{identify.points}{if \code{TRUE}, then identify points interactively.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[base]{par}}).}
  \item{col}{color for points and lines; the default is the \emph{second} entry
    in the current color palette (see \code{\link[base]{palette}}
    and \code{\link[base]{par}}).}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[base]{par}}).}
  \item{lwd}{line width; default is \code{2} (see \code{\link[base]{par}}).}
  \item{main}{title for plot.}
  \item{type}{if \code{"Wang"} use the method of Wang (1985);
    if \code{"Weisberg"} use the method in the Arc software associated with
    Cook and Weisberg (1999).}
  \item{\dots}{arguments to be passed down to \code{av.plot.lm} or
    \code{av.plot.glm}.}
}

\details{
  The function intended for direct use is \code{av.plots} (for which \code{avp}
  is an abbreviation). By default, these functions are used interactively
  through a text menu.
  
  The model can contain factors and interactions. An added-variable plot can be
  drawn for each column of the model matrix, including the constant.
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
  
  Wang, P C. (1985)
  Adding a variable in generalized linear models.
  \emph{Technometrics} \bold{27}, 273--276.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{cr.plots}}, \code{\link{ceres.plots}}, \code{\link{leverage.plots}}}

\examples{
    \dontrun{
av.plots(lm(prestige~income+education+type, data=Duncan))

av.plots(glm(partic != "not.work" ~ hincome + children, 
  data=Womenlf, family=binomial))
    }
}

\keyword{hplot}
\keyword{regression}
