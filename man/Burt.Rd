\name{Burt}
\alias{Burt}
\docType{data}

\title{Fraudulent Data on IQs of Twins Raised Apart}

\description{
The \code{Burt} data frame has 27 rows and 4 columns.
The ``data'' were simply (and notoriously) manufactured.  The
same data are in the dataset ``twins" in the \code{alr3}
package, but with different labels.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{IQbio}{
      IQ of twin raised by biological parents
      }
    \item{IQfoster}{
      IQ of twin raised by foster parents
      }
    \item{class}{
      A factor with levels (note: out of order):
      \code{high};
      \code{low};
      \code{medium}.
      }
    }
  }

\source{
 Burt, C. (1966) 
 The genetic determination of differences in intelligence: 
 A study of monozygotic twins reared together and apart. 
 \emph{British Journal of Psychology} \bold{57}, 137--153.

}


\usage{
Burt
}

\keyword{datasets}
