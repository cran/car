\name{bcPower}
\alias{bcPower}
\alias{bcnPower}
\alias{bcnPowerInverse}
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

bcnPowerInverse(z, lambda, gamma)

yjPower(U, lambda, jacobian.adjusted = FALSE)

basicPower(U,lambda, gamma=NULL)
}
\arguments{
  \item{U}{A vector, matrix or data.frame of values to be transformed}
  \item{lambda}{Power transformation parameter with one element for each 
      column of U, usuallly in the range from \eqn{-2} to \eqn{2}.}
  \item{jacobian.adjusted}{If \code{TRUE}, the transformation is normalized to have
  Jacobian equal to one.  The default \code{FALSE} is almost always appropriate.}
  \item{gamma}{For bcPower or basicPower, the transformation is of U + gamma, where gamma is a positive number called a start that must be large enough so that U + gamma is strictly positive.  For the bcnPower, Box-cox power with negatives allowed, see the details below.}
  \item{z}{a numeric vector the result of a call to \code{bcnPower} with \code{jacobian.adjusted=FALSE}}.
}

\details{
  The Box-Cox 
  family of \emph{scaled power transformations}
  equals \eqn{(x^{\lambda}-1)/\lambda}{(x^(lambda)-1)/lambda}
  for \eqn{\lambda \neq 0}{lambda not equal to 0}, and
  \eqn{\log(x)}{log(x)} if \eqn{\lambda =0}{lambda = 0}.  The \code{bcPower} 
  function computes the scaled power transformation of 
  \eqn{x = U + \gamma}{x = U + gamma}, where \eqn{\gamma}{gamma} 
  is set by the user so \eqn{U+\gamma}{U + gamma} is strictly positive for these
  transformations to make sense.

  The Box-Cox family with negatives allowed was proposed by Hawkins and Weisberg (2017).  It is the Box-Cox power transformation of \deqn{z = .5  (U + \sqrt{U^2 + \gamma^2)})}{z = .5 (U + \sqrt[U^2 + gamma^2])} where for this family \eqn{\gamma}{gamma} is either user selected or is estimated. \code{gamma} must be positive if \eqn{U}{U} includes negative values and non-negative otherwise, ensuring that \eqn{z}{z} is always positive.  The bcnPower transformations behave similarly to the bcPower transformations, and introduce less bias than is introduced by setting the parameter \eqn{\gamma}{gamma} to be non-zero in the Box-Cox family.

The function \code{bcnPowerInverse} computes the inverse of the \code{bcnPower} function, so \code{U = bcnPowerInverse(bcnPower(U, lambda=lam, jacobian.adjusted=FALSE, gamma=gam), lambda=lam, gamma=gam)} is true for any permitted value of \code{gam} and \code{lam}.

  If \code{family="yeo.johnson"} then the Yeo-Johnson transformations are used.
  This is the Box-Cox transformation of \eqn{U+1} for nonnegative values,
  and of \eqn{|U|+1} with parameter \eqn{2-\lambda}{2-lambda} for \eqn{U} negative.

  The basic power transformation returns \eqn{U^{\lambda}}{U^{lambda}} if 
  \eqn{\lambda}{lambda} is not 0, and \eqn{\log(\lambda)}{log(lambda)} 
  otherwise for \eqn{U}{U} strictly positive.

  If \code{jacobian.adjusted} is \code{TRUE}, then the scaled transformations 
  are divided by the
  Jacobian, which is a function of the geometric mean of \eqn{U} for \code{skewPower} and \code{yjPower} and of \eqn{U + gamma} for \code{bcPower}.  With this adjustment, the Jacobian of the transformation is always equal to 1.  Jacobian adjustment facilitates computing the Box-Cox estimates of the transformation parameters.

  Missing values are permitted, and return \code{NA} where ever \code{U} is equal to \code{NA}.

}
\value{
  Returns a vector or matrix of transformed values.
}

\references{
Fox, J. and Weisberg, S. (2019)
\emph{An R Companion to Applied Regression}, Third Edition, Sage.

Hawkins, D. and Weisberg, S. (2017)
Combining the Box-Cox Power and Generalized Log Transformations to Accomodate Nonpositive  Responses In Linear and Mixed-Effects Linear Models \emph{South African Statistics Journal}, 51, 317-328.

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
bcnPower(U, 0, gamma=2)
basicPower(U, lambda = 0, gamma=4)
yjPower(U, 0)
V <- matrix(1:10, ncol=2)
bcPower(V, c(0, 2))
basicPower(V, c(0,1))
}

\keyword{regression}
