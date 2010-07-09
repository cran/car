\name{Greene}
\alias{Greene}
\docType{data}

\title{Refugee Appeals}

\description{
The \code{Greene} data frame has 384 rows and 7 columns.
These are cases filed in 1990, in which refugee claimants rejected
by the Canadian Immigration and Refugee Board asked the Federal
Court of Appeal for leave to appeal the negative ruling of the Board.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{judge}{
      Name of judge hearing case.
      A factor with levels:
      \code{Desjardins}, \code{Heald}, \code{Hugessen}, \code{Iacobucci}, \code{MacGuigan}, 
      \code{Mahoney}, \code{Marceau}, \code{Pratte}, \code{Stone}, \code{Urie}. 
      }
    \item{nation}{
      Nation of origin of claimant.
      A factor with levels:
      \code{Argentina}, \code{Bulgaria}, \code{China}, \code{Czechoslovakia},
      \code{El.Salvador}, \code{Fiji}, \code{Ghana}, \code{Guatemala}, 
      \code{India}, \code{Iran}, \code{Lebanon}, \code{Nicaragua}, 
      \code{Nigeria}, \code{Pakistan}, \code{Poland}, \code{Somalia}, 
      \code{Sri.Lanka}.
      }
    \item{rater}{
      Judgment of independent rater.
      A factor with levels:
      \code{no}, case has no merit;
      \code{yes}, case has some merit (leave to appeal should be granted).
      }
    \item{decision}{
      Judge's decision.
      A factor with levels:
      \code{no}, leave to appeal not granted;
      \code{yes}, leave to appeal granted.
      }
    \item{language}{
      Language of case. 
      A factor with levels:
      \code{English},
      \code{French}.
      }
    \item{location}{
      Location of original refugee claim.
      A factor with levels:
      \code{Montreal},
      \code{other},
      \code{Toronto}.
      }
    \item{success}{
      Logit of success rate, for all cases from the applicant's nation.
      }
    }
  }


\source{
 Personal communication from Ian Greene, Department of
 Political Science, York University.
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
}


\usage{
Greene
}
\keyword{datasets}
