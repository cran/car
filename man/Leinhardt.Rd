\name{Leinhardt}
\alias{Leinhardt}
\docType{data}

\title{Data on Infant-Mortality}

\description{
The \code{Leinhardt} data frame has 105 rows and 4 columns.
The observations are nations of the world around 1970.

}

\format{
  This data frame contains the following columns:
  \describe{
    \item{income}{
      Per-capita income in U. S. dollars.
      }
    \item{infant}{
      Infant-mortality rate per 1000 live births.
      }
    \item{region}{
      A factor with levels:
      \code{Africa}; 
      \code{Americas};
      \code{Asia}, Asia and Oceania;
      \code{Europe}.
      }
    \item{oil}{
      Oil-exporting country. 
      A factor with levels:
      \code{no}, 
      \code{yes}.
      }
    }
  }

\details{
The infant-mortality rate for Jamaica is misprinted in
 Leinhardt and Wasserman; the correct value is given here.
 Some of the values given in Leinhardt and Wasserman do not appear
 in the original New York Times table.
}

\source{
 Leinhardt, S. and Wasserman, S. S. (1979)
 Exploratory data analysis: An introduction to selected methods.
 In Schuessler, K. (Ed.) 
 \emph{Sociological Methodology 1979} Jossey-Bass.
  
 \emph{The New York Times}, 28 September 1975, p. E-3, Table 3.
}

\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}


\usage{
Leinhardt
}
\keyword{datasets}
