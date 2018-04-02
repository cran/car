\name{whichNames}
\alias{which.names}
\alias{whichNames}
\alias{whichNames.data.frame}
\alias{whichNames.default}

\title{Position of Row Names}
\description{
  These functions return the indices of the supplied row names of a data frame
  or names of another object, such as a vector or list. \code{whichNames} is just an alias for \code{which.names}.
}
\usage{
whichNames(names, object, ...)

which.names(names, object, ...)

\method{whichNames}{data.frame}(names, object, ...)

\method{whichNames}{default}(names, object, ...)
}

\arguments{
  \item{names}{a name or character vector of names.}
  \item{object}{a data frame or an object with a names attribute.}
  \item{\dots}{not used.}
}

\value{
Returns the index or indices of \code{names} in row names of the data frame or names of an object of another class. 
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\examples{
whichNames(c('minister', 'conductor'), Duncan)
}
\keyword{utilities}
