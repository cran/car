\name{some}
\alias{some}
\alias{some.data.frame}
\alias{some.matrix}
\alias{some.default}

\title{Sample a Few Elements of an Object}

\description{
 Randomly select a few elements of an object, typically a data frame,
 matrix, vector, or list. If the object is a data frame or a matrix,
 then rows are sampled.
}

\usage{
some(x, ...)

\method{some}{data.frame}(x, n=10, ...)

\method{some}{matrix}(x, n=10, ...)

\method{some}{default}(x, n=10, ...)
}

\arguments{
  \item{x}{the object to be sampled.}
  \item{n}{number of elements to sample.}
  \item{\dots}{arguments passed down.}
}

\value{
    Sampled elements or rows. 
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}

\note{These functions are adapted from \code{head} and \code{tail} 
in the \code{utils} package.
}


\seealso{\code{\link[utils]{head}}, \code{\link[utils]{tail}}.}

\examples{
some(Duncan)
}
\keyword{utilities}
