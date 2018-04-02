\name{sigmaHat}
\alias{sigmaHat}
\alias{sigmaHat.default}
\alias{sigmaHat.glm}
\alias{sigmaHat.lm}

\title{ Return the scale estimate for a regression model }
\description{
  This function provides a consistent method to return the estimated scale
  from a linear, generalized linear, nonlinear, or other model.
}
\usage{
sigmaHat(object)
}
\arguments{
  \item{object}{A regression object of type lm, glm or nls}
}
\details{
For an lm or nls object, the returned quantity is the square root of the
estimate of \eqn{\sigma^2}{sigma^2}.  For a glm object, the returned quantity is
the square root of the estimated dispersion parameter.
}
\value{
  A nonnegative number
}
\author{Sanford Weisberg, \email{sandy@umn.edu} }

\examples{
m1 <- lm(prestige ~ income + education, data=Duncan)
sigmaHat(m1)
}
\keyword{ regression }% at least one, from doc/KEYWORDS
