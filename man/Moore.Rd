\name{Moore}
\alias{Moore}
\docType{data}

\title{Status, Authoritarianism, and Conformity}

\description{
The \code{Moore} data frame has 45 rows and 4 columns.
The data are for subjects in a social-psychological experiment,
who were faced with manipulated disagreement from a partner of either
of low or high status. The subjects could either conform to the
partner's judgment or stick with their own judgment.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{partner.status}{
      Partner's status.  A factor with levels:
      \code{high},
      \code{low}.
      }
    \item{conformity}{
      Number of conforming responses in 40 critical trials.
      }
    \item{fcategory}{
      F-Scale Categorized.
      A factor with levels (note levels out of order):
      \code{high},
      \code{low},
      \code{medium}.
      }
    \item{fscore}{
      Authoritarianism: F-Scale score.
      }
    }
  }
 
\source{
 Moore, J. C., Jr. and Krupat, E. (1971) 
 Relationship between source status, authoritarianism and conformity in a
 social setting. \emph{Sociometry} \bold{34}, 122--134.  
 
 Personal communication
 from J. Moore, Department of Sociology, York University.
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}


\usage{
Moore
}
\keyword{datasets}
