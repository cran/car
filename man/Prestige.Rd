\name{Prestige}
\alias{Prestige}
\docType{data}

\title{Prestige of Canadian Occupations}

\description{
The \code{Prestige} data frame has 102 rows and 6 columns.
The observations are occupations.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{education}{
      Average education of occupational incumbents, years, in 1971. 
      }
    \item{income}{
      Average income of incumbents, dollars, in 1971.
      }
    \item{women}{
      Percentage of incumbents who are women.
      }
    \item{prestige}{
      Pineo-Porter prestige score for occupation, from a social survey
      conducted in the mid-1960s.
      }
    \item{census}{
      Canadian Census occupational code.
      }
    \item{type}{
      Type of occupation.
      A factor with levels (note: out of order):
      \code{bc}, Blue Collar;
      \code{prof}, Professional, Managerial, and Technical;
      \code{wc}, White Collar. 
      }
    }
  }

\source{
 Canada (1971)
 \emph{Census of Canada}. Vol. 3, Part 6. Statistics Canada [pp. 19-1--19-21].
 
 Personal communication from B. Blishen, W. Carroll, and
 C. Moore, Departments of Sociology, York University and
 University of Victoria.
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}


\usage{
Prestige
}

\keyword{datasets}
