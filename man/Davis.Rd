\name{Davis}
\alias{Davis}
\non_function{}

\title{Self-Reports of Height and Weight}

\description{
The \code{Davis} data frame has 200 rows and 5 columns.
The subjects were men and women engaged in regular exercise.
There are some missing data.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{sex}{
      A factor with levels:
      \code{F}, male;
      \code{M}, female.
      }
    \item{weight}{
      Measured weight in kg. 
      }
    \item{height}{
      Measured height in cm. 
      }
    \item{repwt}{
      Reported weight in kg. 
      }
    \item{repht}{
      Reported height in cm. 
      }
    }
  }


\source{
 Personal communication from C. Davis, Departments of
 Physical Education and Psychology, York University.
}

\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.
}


\usage{
data(Davis)
}

\keyword{datasets}
