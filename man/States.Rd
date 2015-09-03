\name{States}
\alias{States}
\docType{data}

\title{Education and Related Statistics for the U.S. States}

\description{
  The \code{States} data frame has 51 rows and 8 columns.
  The observations are the U. S. states and Washington, D. C.
}

\format{
  This data frame contains the following columns:
  \describe{
    \item{region}{
      U. S. Census regions.
      A factor with levels:
      \code{ENC}, East North Central;
      \code{ESC}, East South Central;
      \code{MA}, Mid-Atlantic;
      \code{MTN}, Mountain; 
      \code{NE}, New England;
      \code{PAC}, Pacific; 
      \code{SA}, South Atlantic;
      \code{WNC}, West North Central; 
      \code{WSC}, West South Central. 
      }
    \item{pop}{
      Population: in 1,000s.
      }
    \item{SATV}{
      Average score of graduating high-school 
      students in the state on the \emph{verbal} component of the 
      Scholastic Aptitude Test 
      (a standard university admission exam).
      }
    \item{SATM}{
      Average score of graduating high-school 
      students in the state on the \emph{math} component of the 
      Scholastic Aptitude Test.
      }
    \item{percent}{
       Percentage of graduating high-school 
      students in the state who took the SAT exam.
      }
    \item{dollars}{
       State spending on public education, 
      in \$1000s per student.
      }
    \item{pay}{
       Average teacher's salary in the state, 
      in $1000s.
      }
    }
  }

\source{
 United States (1992)
 \emph{Statistical Abstract of the United States.} Bureau of the Census.
}

\references{
 Moore, D. (1995)
 \emph{The Basic Practice of Statistics}. Freeman, Table 2.1.
}

\usage{
States
}

\keyword{datasets}
