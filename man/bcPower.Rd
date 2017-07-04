\name{bcPower}
\alias{bcPower}
\alias{bcnPower}
\alias{yjPower}
\alias{basicPower}

\title{Box-Cox, Box-Cox with Negatives Allowed, Yeo-Johnson and Basic Power Transformations}
\description{
Transform the elements of a vector or columns of a matrix using, the Box-Cox, Box-Cox with negatives allowed,
Yeo-Johnson, or simple power transformations.
}
\usage{
bcPower(U, lambda, jacobian.adjusted=FALSE, gamma=NULL)

bcnPower(U, lambda, jacobian.adjusted = FALSE, gamma)

yjPower(U, lambda, jacobian.adjusted = FALSE)

basicPower(U,lambda, gamma=NULL)
}
\arguments{
  \item{U}{A vector, matrix or data.frame of values to be transformed}
  \item{lambda}{Power transformation parameter with one element for eacul column of U, usuallly in
  the range from \eqn{-2} to \eqn{2}, or if \code{U}}
  \item{jacobian.adjusted}{If \code{TRUE}, the transformation is normalized to have
  Jacobian equal to one.  The default \code{FALSE} is almost always appropriate}
  \item{gamma}{For bcPower or basicPower, the transformation is of U + gamma, where gamma is a positive number called a start that must be large enough so that U + gamma is strictly positive.  For the bcnPower, Box-cox power with negatives allowed, see the details below.}
}

\details{
  The Box-Cox family of \emph{scaled power transformations}
  equals \eqn{((U + \gamma)^{\lambda}-1)/\lambda}{(U^(lambda)-1)/lambda}
  for \eqn{\lambda \neq 0}{lambda not equal to zero}, and
  \eqn{\log(U)}{log(U)} if \eqn{\lambda =0}{lambda = 0}.  If \eqn{\gamma}{gamma} is not specified, it is set equal to zero.  \code{U + gamma}
  must be strictly positive to use this family.

  The Box-Cox family with negatives allowed was proposed by Hawkins and Weisberg (2017).  It is the Box-Cox power transformation of \eqn{z = .5 * (y + (y^2 + \gamma^2)^{1/2})}, where \eqn{\gamma}{gamma} is strictly positive if \eqn{y}{y} includes negative values and non-negative otherwise.  The value of \eqn{z}{z} is always positive.  The bcnPower transformations behave very similarly to the bcPower transformations, including much less bias than is introduced by setting the parameter \eqn{\gamma}{gamma} to be non-zero in the Box-Cox family.

  If \code{family="yeo.johnson"} then the Yeo-Johnson transformations are used.
  This is the Box-Cox transformation of \eqn{U+1} for nonnegative values,
  and of \eqn{|U|+1} with parameter \eqn{2-\lambda}{2-lambda} for \eqn{U} negative.

  The basic power transformation returns \eqn{U^{\lambda}} if \eqn{\lambda}
  is not zero, and \eqn{\log(\lambda)} otherwise for \eqn{U}{U} strictly positive.

  If \code{jacobian.adjusted} is \code{TRUE}, then the scaled transformations are divided by the
  Jacobian, which is a function of the geometric mean of \eqn{U} for \code{skewPower} and \code{yjPower} and of \eqn{U + gamma} for \code{bcPower}.  With this adjustment, the Jacobian of the transformation is always equal to 1.

  Missing values are permitted, and return \code{NA} where ever \code{U} is equal to \code{NA}.

}
\value{
  Returns a vector or matrix of transformed values.
}

\references{
Fox, J. and Weisberg, S. (2011)
\emph{An R Companion to Applied Regression}, Second Edition, Sage.

Hawkins, D. and Weisberg, S. (2017)
Combining the Box-Cox Power and Generalized Log Transformations to Accomodate Negative Responses
In Linear and Mixed-Effects Linear Models,
submitted for publication.

Weisberg, S. (2014)  \emph{Applied Linear Regression}, Fourth Edition, Wiley
Wiley, Chapter 7.

Yeo, In-Kwon and Johnson, Richard (2000) A new family of power
transformations to improve normality or symmetry.  \emph{Biometrika}, 87,
954-959.
}
\author{ Sanford Weisberg, <sandy@umn.edu> }

\seealso{\code{\link{powerTransform}}, \code{\link{testTransform}}}
\examples{
U <- c(NA, (-3:3))
\dontrun{bcPower(U, 0)}  # produces an error as U has negative values
bcPower(U, 0, gamma=4)
bcPower(U, .5, jacobian.adjusted=TRUE, gamma=4)
basicPower(U, lambda = 0, gamma=4)
yjPower(U, 0)
V <- matrix(1:10, ncol=2)
bcPower(V, c(0, 2))
basicPower(V, c(0,1))
}

\keyword{regression}
