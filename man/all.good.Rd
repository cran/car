\name{all.good}
\alias{all.good}

\title{Function to Test for Non-Missing Observations}
\description{
 This function finds observations that are non-missing for 
 \emph{all} of its arguments.
}
\usage{
all.good(...)
}

\arguments{
  \item{\dots}{one or more vectors or factors of the same length; may
  include matrices with the same number of rows.}
}
\
\value{
  A logical vector with TRUE for observations that are non-missing
  for \emph{all} arguments, normally used to index valid observations.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
w <- matrix(1:20, 10, 2)
x <- 1:10
y <- factor(rep(c('a', 'b'), 5))
z <- letters[1:10]
w[1,1] <- NA
x[3] <- NA
y[5] <- NA
z[7] <- NA
all.good(w, x, y, z)
##  [1] FALSE  TRUE FALSE  TRUE FALSE  TRUE FALSE  TRUE  TRUE  TRUE
}

\keyword{utilities}
