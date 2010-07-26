\name{which.names}
\alias{which.names}
\alias{whichNames}

\title{Position of Row Names}
\description{
  These functions return the indices of row names in a data frame
  or a vector of names. \code{whichNames} is just an alias for \code{which.names}.
}
\usage{
which.names(names, object)
whichNames(...)
}

\arguments{
  \item{names}{a name or character vector of names.}
  \item{object}{a data frame or character vector of (row) names.}
  \item{\dots}{arguments to be passed to \code{which.names}.}
}

\value{
Returns the index or indices of \code{names} within \code{object}.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}

\examples{
which.names(c('minister', 'conductor'), Duncan)
## [1]  6 16
}
\keyword{utilities}
