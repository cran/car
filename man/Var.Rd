\name{Var}
\alias{Var}
\alias{Var.default}
\alias{Var.lm}
\alias{Var.glm}

\title{Variance-Covariance Matrices (deprecated)}
\description{
  Computes variance-covariance matrices or variances for model objects or data.
  The default method uses the function \code{var}. 
  
  These functions are now deprecated; instead, use
  the \code{vcov} function, now in the base package. Note that \code{vcov} has no \code{diagonal}
  argument and no default method.
}

\usage{
Var(object, ...)

\method{Var}{default}(object, diagonal=FALSE, ...)

\method{Var}{lm}(object, diagonal=FALSE, ...)

\method{Var}{glm}(object, diagonal=FALSE, ...)
}

\arguments{
  \item{object}{an object for which the covariance matrix is desired.}
  \item{\dots}{arguments to be passed to \code{var} (e.g., \code{na.rm}).}
  \item{diagonal}{if \code{TRUE}, return only the variances.}
}

\value{
  A variance-covariance matrix or a vector of variances.
  }

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[stats]{var}}}

\examples{
data(Davis)
attach(Davis)
Var(cbind(weight, repwt), na.rm=TRUE)
##          weight    repwt
## weight 233.8781 176.1014
## repwt  176.1014 189.7966

Var(lm(weight~repwt))
##              (Intercept)        repwt
##  (Intercept)   9.2228211 -0.134640952
##  repwt        -0.1346410  0.002051736
}

\keyword{misc}
