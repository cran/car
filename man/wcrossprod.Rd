\name{wcrossprod}
\alias{wcrossprod}
\title{
Weighted Matrix Crossproduct
}
\description{
Given matrices \code{x} and \code{y} as arguments and an optional matrix or vector of weights, \code{w}, 
return a weighted matrix cross-product, \code{t(x) w y}.
If no weights are supplied, or the weights are constant, the function uses
\code{\link[base]{crossprod}} for speed.

}
\usage{
wcrossprod(x, y, w)
}

\arguments{
  \item{x,y}{
x, y numeric matrices; \code{missing(y)} is taken to be the same matrix as x. 
Vectors are promoted to single-column or single-row matrices, depending on the context. 
}
  \item{w}{
A numeric vector or matrix of weights, conformable with \code{x} and \code{y}.
}
}

\value{
A numeric matrix, with appropriate dimnames taken from \code{x} and \code{y}.

}

\author{
Michael Friendly,
John Fox \email{jfox@mcmaster.ca}
}


\seealso{
\code{\link[base]{crossprod}}
}
\examples{
set.seed(12345)
n <- 24
drop <- 4
sex <- sample(c("M", "F"), n, replace=TRUE)
x1 <- 1:n
x2 <- sample(1:n)
extra <- c( rep(0, n - drop), floor(15 + 10 * rnorm(drop)) )
y1 <- x1 + 3*x2 + 6*(sex=="M") + floor(10 * rnorm(n)) + extra
y2 <- x1 - 2*x2 - 8*(sex=="M") + floor(10 * rnorm(n)) + extra
# assign non-zero weights to 'dropped' obs
wt <- c(rep(1, n-drop), rep(.2,drop))

X <- cbind(x1, x2)
Y <- cbind(y1, y2)
wcrossprod(X)
wcrossprod(X, w=wt)

wcrossprod(X, Y)
wcrossprod(X, Y, w=wt)

wcrossprod(x1, y1)
wcrossprod(x1, y1, w=wt)

}

\keyword{ array}
\keyword{algebra}
