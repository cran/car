\name{influence.plot}
\alias{influence.plot}

\title{Regression Influence Plot}
\description{
  This function creates a "bubble" plot of studentized residuals by hat values, with the areas of the
  circles representing the observations proportional to Cook's distances. Vertical reference lines are drawn
  at twice and three times the average hat value, horizontal reference lines at -2, 0, and 2 on the
  studentized-residual scale.
}
\usage{
influence.plot(model, scale = 10, col = c(1, 2), labels = names(rstud), ...)
}

\arguments{
  \item{model}{a linear or generalized-linear model.}
  \item{scale}{a factor to adjust the size of the circles.}
  \item{col}{colors for plotting points that do not and do have noteworthy Cook's distances.}
  \item{labels}{if \code{FALSE} do not identify points interactively with the mouse; otherwise
    a vector of observation labels.}
  \item{\dots}{arguments to pass to the \code{plot} function.}
}

}
\value{
  Returns the indices of identified points.
}

\references{
    J. Fox (2002) 
    \emph{An R and S-PLUS Companion to Applied Regression}. Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}


\seealso{\code{\link[car]{cookd}}, \code{\link[base]{rstudent}}, \code{\link[base]{hatvalues}}}

\examples{
    \dontrun{
    data(Duncan)
    influence.plot(lm(prestige ~ income + education, data=Duncan))
    }
}

\keyword{regression}
