\name{Var}
\alias{Var}
\alias{Var.default}
\alias{Var.lm}
\alias{Var.glm}

\title{Variance-Covariance Matrices}
\description{
  Computes variance-covariance matrices or variances for model objects or data.
  The default methods uses the function \code{var}.
}

\usage{
Var(object, ...)

Var.default(object, diagonal=F, ...)

Var.lm(object, diagonal=F)

Var.glm(object, diagonal=F)
}

\arguments{
  \item{object}{an object for which the covariance matrix is desired.}
  \item{\dots}{arguments to be passed to \code{var} (e.g., \code{na.rm}).}
  \item{diagonal}{if \code{TRUE}, return only the variances.}
}

\value{
  A variance-covariance matrix or a vector of variances.

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[base]{var}}}

\examples{
data(Davis)
attach(Davis)
Var(cbind(weight, repwt), na.rm=T)
##          weight    repwt
## weight 233.8781 176.1014
## repwt  176.1014 189.7966

Var(lm(weight~repwt))
##              (Intercept)        repwt
##  (Intercept)   9.2228211 -0.134640952
##  repwt        -0.1346410  0.002051736
}

\keyword{misc}
