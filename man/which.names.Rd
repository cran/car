\name{which.names}
\alias{which.names}

\title{Position of Row Names}
\description{
  This function returns the indices of row names in a data frame
  or a vector of names.
}
\usage{
which.names(names, object)
}

\arguments{
  \item{names}{a name or character vector of names.}
  \item{object}{a data frame or character vector of (row) names.}
}

\value{
Returns the index or indices of \code{names} within \code{object}.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
which.names(c('minister', 'conductor'), Duncan)
## [1]  6 16
}
\keyword{utilities}
