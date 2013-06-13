\name{Vocab}
\alias{Vocab}
\docType{data}

\title{Vocabulary and Education}

\description{
The \code{Vocab} data frame has 21,638 rows and 5 columns.
The observations are respondents to U.S. General Social Surveys, 1972-2004.
}

\format{
  This data frame contains the following columns:
  \describe{
  	\item{year}{Year of the survey.}
  	\item{sex}{Sex of the respondent, \code{Female} or \code{Male}.}
    \item{education}{
      Education, in years.
      }
    \item{vocabulary}{
      Vocabulary test score: number correct on a 10-word test.
      }
    }
  }

\source{
 National Opinion Research Center
 \emph{General Social Survey.} 
 GSS Cumulative Datafile 1972-2004, downloaded from \url{http://sda.berkeley.edu/archive.htm}.
}

\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}


\usage{
Vocab
}

\keyword{datasets}
