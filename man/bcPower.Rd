\name{bcPower}
\alias{bcPower}
\alias{yjPower}
\alias{basicPower}

\title{Box-Cox and Yeo-Johnson Power Transformations}
\description{
Transform the elements of a vector using, the Box-Cox,
Yeo-Johnson, or simple power transformations.
}
\usage{
bcPower(U, lambda, jacobian.adjusted = FALSE)

yjPower(U, lambda, jacobian.adjusted = FALSE)

basicPower(U,lambda)
}
\arguments{
  \item{U}{A vector, matrix or data.frame of values to be transformed}
  \item{lambda}{The one-dimensional transformation parameter, usually in 
  the range from \eqn{-2} to \eqn{2}, or if \code{U} is a matrix or data frame, a vector of length
  \code{ncol(U)} of transformation parameters}
  \item{jacobian.adjusted}{If \code{TRUE}, the transformation is normalized to have
  Jacobian equal to one.  The default is \code{FALSE}.}
}

\details{
  The Box-Cox family of \emph{scaled power transformations} 
  equals \eqn{(U^{\lambda}-1)/\lambda}{(U^(lambda)-1)/lambda} 
  for \eqn{\lambda \neq 0}{lambda not equal to zero}, and 
  \eqn{\log(U)}{log(U)} if \eqn{\lambda =0}{lambda = 0}.
  
  If \code{family="yeo.johnson"} then the Yeo-Johnson transformations are used.  
  This is the Box-Cox transformation of \eqn{U+1} for nonnegative values, 
  and of \eqn{|U|+1} with parameter \eqn{2-\lambda}{2-lambda} for \eqn{U} negative.

  If \code{jacobian.adjusted} is \code{TRUE}, then the scaled transformations are divided by the 
  Jacobian, which is a function of the geometric mean of \eqn{U}.
  
  The basic power transformation returns \eqn{U^{\lambda}} if \eqn{\lambda}
  is not zero, and \eqn{\log(\lambda)} otherwise.
  
  Missing values are permitted, and return \code{NA} where ever \code{U}is equal to \code{NA}.
  
}
\value{
  Returns a vector or matrix of transformed values.
}

\references{
Fox, J. and Weisberg, S. (2011) 
\emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
Weisberg, S. (2005)  \emph{Applied Linear Regression}, Third Edition.
Wiley, Chapter 7.

Yeo, In-Kwon and Johnson, Richard (2000) A new family of power
transformations to improve normality or symmetry.  \emph{Biometrika}, 87,
954-959.
}
\author{ Sanford Weisberg, <sandy@stat.umn.edu> }

\seealso{\code{\link{powerTransform}} }
\examples{
U <- c(NA, (-3:3))
\dontrun{bcPower(U, 0)}  # produces an error as U has negative values
bcPower(U+4,0)
bcPower(U+4, .5, jacobian.adjusted=TRUE)
yjPower(U, 0)
yjPower(U+3, .5, jacobian.adjusted=TRUE)
V <- matrix(1:10, ncol=2)
bcPower(V, c(0,1))
#basicPower(V, c(0,1))
}

\keyword{regression}
