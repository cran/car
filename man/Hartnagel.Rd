\name{Hartnagel}
\alias{Hartnagel}
\docType{data}

\title{Canadian Crime-Rates Time Series}

\description{
The \code{Hartnagel} data frame has 38 rows and 7 columns.
The data are an annual time-series from 1931 to 1968. There are
some missing data.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{year}{
      1931--1968.
      }
    \item{tfr}{
      Total fertility rate per 1000 women.
      }
    \item{partic}{
      Women's labor-force participation rate per 1000.
      }
    \item{degrees}{
      Women's post-secondary degree rate per 10,000.
      }
    \item{fconvict}{
      Female indictable-offense conviction rate per 100,000.
      }
    \item{ftheft}{
      Female theft conviction rate per 100,000.
      }
    \item{mconvict}{
      Male indictable-offense conviction rate per 100,000.
      }
    \item{mtheft}{
      Male theft conviction rate per 100,000.
      }
    }
  }

\details{
 The post-1948 crime rates have been adjusted to account for
 a difference in method of recording. Some of your results will differ
 in the last decimal place from those in Table 14.1 of Fox (1997) due
 to rounding of the data. Missing values for 1950 were interpolated.
}


\source{
 Personal communication from T. Hartnagel, 
 Department of Sociology, University of Alberta.
}

\references{
  Fox, J., and Hartnagel, T. F (1979)
  Changing social roles and female crime in Canada:
  A time series analysis.
  \emph{Canadian Review of Sociology and Anthroplogy}, 
  \bold{16}, 96--104.
  
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
}

\usage{
Hartnagel
}

\keyword{datasets}
