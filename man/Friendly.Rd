\name{Friendly}
\alias{Friendly}
\docType{data}

\title{Format Effects on Recall}

\description{
 The \code{Friendly} data frame has 30 rows and 2 columns.
 The data are from an experiment on subjects' ability to remember words
 based on the presentation format.

}

\format{
  This data frame contains the following columns:
    \describe{
      \item{condition}{
        A factor with levels:
          \code{Before}, Recalled words presented before others;
          \code{Meshed}, Recalled words meshed with others; 
          \code{SFR}, Standard free recall.
          }
      \item{correct}{
        Number of words correctly recalled, out of 40 on final trial of the experiment.
        }
     }
  }


\source{
 Friendly, M. and Franklin, P. (1980) 
 Interactive presentation in multitrial free recall. 
 \emph{Memory and Cognition}
 \bold{8} 265--270 [Personal communication from M. Friendly].
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}


\usage{
Friendly
}

\keyword{datasets}
