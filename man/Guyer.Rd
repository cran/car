\name{Guyer}
\alias{Guyer}
\docType{data}

\title{Anonymity and Cooperation}

\description{
The \code{Guyer} data frame has 20 rows and 3 columns.
The data are from an experiment in which four-person groups
played a prisoner's dilemma game for 30 trails, each person
making either a cooperative or competitive choice on each
trial. Choices were made either anonymously or in public;
groups were composed either of females or of males.
The observations are 20 groups.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{cooperation}{
      Number of cooperative choices (out of 120 in all).
      }
    \item{condition}{
      A factor with levels:
      \code{A}, Anonymous;
      \code{P}, Public-Choice.
      }
    \item{sex}{
      Sex.  A factor with levels:
      \code{F}, Female;
      \code{M}, Male.
      }
    }
  }
 
\source{
 Fox, J. and Guyer, M. (1978)
 Public choice and cooperation in n-person prisoner's dilemma.
 \emph{Journal of Conflict Resolution} \bold{22}, 469--481.
}

\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}

\usage{
Guyer
}

\keyword{datasets}
