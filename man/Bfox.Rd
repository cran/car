\name{Bfox}
\alias{Bfox}
\docType{data}

\title{Canadian Women's Labour-Force Participation}
\description{
The \code{Bfox} data frame has 30 rows and 7 columns.
Time-series data on Canadian women's labor-force participation,
1946--1975.

}

\format{
  This data frame contains the following columns:
  \describe{
    \item{partic}{
      Percent of adult women in the workforce.
      }
    \item{tfr}{
      Total fertility rate: expected births to a cohort of 1000
      women at current age-specific fertility rates.
      }
    \item{menwage}{
      Men's average weekly wages, in constant 1935 dollars and
      adjusted for current tax rates.
      }
    \item{womwage}{
      Women's average weekly wages.
      }
    \item{debt}{
      Per-capita consumer debt, in constant dollars.
      }
    \item{parttime}{
      Percent of the active workforce working 34 hours per week or
      less.
      }
    }
  }
  
\section{Warning}{
The value of \code{tfr} for 1973 is misrecorded as 2931; it should be 1931.
}


\source{
   Fox, B. (1980) 
   \emph{Women's Domestic Labour and their Involvement 
   in Wage Work.} Unpublished doctoral dissertation, p. 449.
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
}

\usage{
Bfox
}

\keyword{datasets}
