\name{Duncan}
\alias{Duncan}
\docType{data}

\title{Duncan's Occupational Prestige Data}

\description{
  The \code{Duncan} data frame has 45 rows and 4 columns.
  Data on the prestige and other characteristics of 45 U. S. occupations in 1950.
  }
  
\format{
  This data frame contains the following columns:
  \describe{
    \item{type}{ Type of occupation.
      A factor with the following levels:
        \code{prof}, professional and managerial;
        \code{wc}, white-collar;
        \code{bc}, blue-collar.
      }
    \item{income}{
      Percent of males in occupation earning $3500 or more in 1950.
      }
    \item{education}{
      Percent of males in occupation in 1950 who were high-school graduates.
      }
    \item{prestige}{
      Percent of raters in NORC study rating occupation as excellent or good in prestige.
      }
    }
  }

\source{
  Duncan, O. D. (1961) 
  A socioeconomic index for all occupations.
  In Reiss, A. J., Jr. (Ed.)
  \emph{Occupations and Social Status.} Free Press 
  [Table VI-1].
  }
 
\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}

\usage{
Duncan
}

\keyword{datasets}
