\name{influencPlot}
\alias{influencePlot}
\alias{influencePlot.lm}
\alias{influence.plot}

\title{Regression Influence Plot}
\description{
  This function creates a "bubble" plot of studentized residuals by hat values, with the areas of the
  circles representing the observations proportional to Cook's distances. Vertical reference lines are drawn
  at twice and three times the average hat value, horizontal reference lines at -2, 0, and 2 on the
  studentized-residual scale.
}

\note{
  This function used to be named \code{influence.plot}; the name was changed to avoid confusion with
  the \code{influence} generic function.
}

\usage{
influencePlot(model, ...)

\method{influencePlot}{lm}(model, scale=10, col=c(1,2), identify=c(TRUE, FALSE, "auto"),
		labels=names(rstud), cex.identify=par("cex"), col.identify=par("col"), ...)
}

\arguments{
  \item{model}{a linear or generalized-linear model.}
  \item{scale}{a factor to adjust the size of the circles.}
  \item{col}{colors for plotting points that do not and do have noteworthy Cook's distances.}
  \item{identify}{identify points; if \code{TRUE}, the default, identify points interactively;
  	if \code{"auto"} then points with large Cook's distances will automatically be identified.}
  \item{labels}{a vector of observation labels.}
  \item{cex.identify, col.identify}{for point labels.}
  \item{\dots}{arguments to pass to the \code{plot} function.}
}

\value{
  Returns the indices of identified points.
}

\references{
    J. Fox (2002) 
    \emph{An R and S-PLUS Companion to Applied Regression}. Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}


\seealso{\code{\link[car]{cookd}}, \code{\link[stats]{rstudent}}, \code{\link[stats]{hatvalues}}}

\examples{
    influencePlot(lm(prestige ~ income + education, data=Duncan), 
    	identify="auto")
}

\keyword{regression}
