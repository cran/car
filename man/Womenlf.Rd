\name{Womenlf}
\alias{Womenlf}
\docType{data}

\title{Canadian Women's Labour-Force Participation}

\description{
  The \code{Womenlf} data frame has 263 rows and 4 columns.
  The data are from a 1977 survey of the Canadian population.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{partic}{
      Labour-Force Participation.
      A factor with levels (note: out of order):
      \code{fulltime}, Working full-time;
      \code{not.work}, Not working outside the home;
      \code{parttime}, Working part-time.
      }
    \item{hincome}{
      Husband's income, $1000s.
      }
    \item{children}{
      Presence of children in the household.
      A factor with levels:
      \code{absent}, 
      \code{present}.
      }
    \item{region}{
      A factor with levels:
      \code{Atlantic}, Atlantic Canada;
      \code{BC}, British Columbia;
      \code{Ontario}; 
      \code{Prairie}, Prairie provinces;
      \code{Quebec}.
      }
    }
  }

\source{
 \emph{Social Change in Canada Project.} York Institute for
 Social Research.
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}

\usage{
Womenlf
}

\keyword{datasets}
