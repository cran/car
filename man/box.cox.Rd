\name{box.cox}
\alias{box.cox}
\alias{bc}

\title{Box-Cox Family of Transformations}
\description{
  Compute the Box-Cox power transformation of a variable.
}
\usage{
box.cox(x, p, start=0)

bc(x, p, ...)
}

\arguments{
  \item{x}{numeric vector to transform.}
  \item{p}{power (0 = log); if \code{p} is a vector then a matrix
    of transformed values with columns labelled by powers will be
    returned.}
  \item{start}{constant to be added to each value of \code{x} prior
    to transformation.}
  \item{...}{argument passed down.}
}
\details{
  Computes \eqn{x^{\prime }=(x^p-1)/p}{x' = (x^p - 1)/p} 
  for \eqn{p \neq 0}{p != 0} and \eqn{x^{\prime }=\log x}{x' = log(x)} 
  for \eqn{p=0}{p = 0}.
  
  The values of \code{x} must all be positive; if not, a \code{start} should
  be added to each value to make all the values positive. The function will
  automatically compute the \code{start} and print a warning, if necessary. 
  
  \code{bc} is just an abbreviation for \code{box.cox}.
}
\value{
  a vector or matrix of transformed values.
}
\references{
  Atkinson, A. C. (1985)
  \emph{Plots, Transformations, and Regression} Oxford.
  
  Box, G. E. P. and Cox, D. R. (1964)
  An analysis of transformations.
  \emph{JRSS B} \bold{26}, 211--246.
  
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage. 
}

\author{John Fox \email{jfox@mcmaster.ca}}

\section{Warning}{These functions \emph{do not} compute the maximum-likelihood estimate for a Box-Cox
normalizing transformation. See \code{\link{box.cox.powers}} for estimating unconditional
univariate and multivariate Box-Cox transformations, and \code{boxcox} in the
\code{MASS} package for estimating the Box-Cox transformation of the response in a
linear model.}

\seealso{\code{\link[MASS]{boxcox}}, \code{\link{box.cox.var}}, 
  \code{\link{box.cox.powers}}, \code{\link{box.cox.axis}}}

\examples{
box.cox(1:10, 2)
## [1]  0.0  1.5  4.0  7.5 12.0 17.5 24.0 31.5 40.0 49.5

box.cox(1:5, c(0,2))
##              0    2
## [1,] 0.0000000  0.0
## [2,] 0.6931472  1.5
## [3,] 1.0986123  4.0
## [4,] 1.3862944  7.5
## [5,] 1.6094379 12.0

box.cox(-5:5, 2)
##  [1]  0.0  1.5  4.0  7.5 12.0 17.5 24.0 31.5 40.0 49.5 60.0
## Warning message: 
## start =  6 added to data prior to transformation in: box.cox(-5:5, 2)

options(digits=4)
box.cox(-5:5, 0, start=6)
##  [1] 0.0000 0.6931 1.0986 1.3863 1.6094 1.7918 1.9459 2.0794 2.1972
## [10] 2.3026 2.3979
}

\keyword{manip}
\keyword{regression}
