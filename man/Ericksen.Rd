\name{Ericksen}
\alias{Ericksen}
\docType{data}

\title{The 1980 U.S. Census Undercount}

\description{
The \code{Ericksen} data frame has 66 rows and 9 columns.
The observations are 16 large cities, the remaining parts of
the states in which these cities are located, and the other
U. S. states.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{minority}{
      Percentage black or Hispanic.
      }
    \item{crime}{
      Rate of serious crimes per 1000 population.
      }
    \item{poverty}{
      Percentage poor.
      }
    \item{language}{
      Percentage having difficulty speaking or writing English.
      }
    \item{highschool}{
      Percentage age 25 or older who had not finished highschool.
      }
    \item{housing}{
      Percentage of housing in small, multiunit buildings.
      }
    \item{city}{A factor with levels: 
      \code{city}, major city; 
      \code{state}, state or state-remainder.
      }
    \item{conventional}{
      Percentage of households counted by conventional personal enumeration.
      }
    \item{undercount}{
      Preliminary estimate of percentage undercount.
      }
    }
  }

\source{
Ericksen, E. P., Kadane, J. B. and Tukey, J. W. (1989)
 Adjusting the 1980 Census of Population and Housing. 
 \emph{Journal of the American Statistical Association} \bold{84}, 927--944 [Tables 7 and 8].
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}

\usage{
Ericksen
}

\keyword{datasets}
