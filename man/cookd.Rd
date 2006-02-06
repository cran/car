\name{Cook's Distances}

\alias{cookd}

\title{Cook's Distances for Linear and Generalized Linear Models}
\description{
  This function now simply calls \code{cooks.distance} in the \code{base} package.
}

\usage{
cookd(model, ...)
}

\arguments{
  \item{model}{\code{lm} or \code{glm} model object.}
  \item{...}{other arguments to be passed to \code{cooks.distance}.}
}

\details{
  Cook's distances for generalized linear
  models are approximations, as described in Williams (1987) (except that the Cook's distances are
  scaled as \emph{F} rather than as chi-square values).
  
  This function is retained primarily for consistency with \emph{An R and S-PLUS Companion
  to Applied Regression.} Other deletion diagnostics formerly in the \code{car} package have
  been rewritten and moved to the \code{base} package; these include \code{influence}, 
  \code{rstudent}, \code{hatvalues}, \code{dfbeta}, and \code{dfbetas}.
}

\value{
 \code{cookd} returns a vector with one entry for each observation.
}

\references{
  Cook, R. D. and Weisberg, S. (1984)
  \emph{Residuals and Influence in Regression.} Wiley.
  
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
  
  Williams, D. A. (1987)
  Generalized linear model diagnostics using the deviance and single
  case deletions. \emph{Applied Statistics} \bold{36}, 181--191.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[base:cooks.distance]{cooks.distance}}, \code{\link[base:influence]{influence}}}

\examples{
plot(cookd(lm(prestige ~ income + education, data=Duncan)))
}

\keyword{regression}
\keyword{models}
