\name{Salaries}
\Rdversion{1.1}
\alias{Salaries}
\docType{data}
\title{
Salaries for Professors
}

\description{
The 2008-09 nine-month academic
salary for Assistant Professors, Associate Professors and Professors in a
college in the U.S. The data were collected as part of the
on-going effort of the college's administration to monitor salary differences
between male and female faculty members. 
}

\usage{Salaries}
\format{
  A data frame with 397 observations on the following 6 variables.
  \describe{
    \item{\code{rank}}{a factor with levels \code{AssocProf} \code{AsstProf} \code{Prof}}
    \item{\code{discipline}}{a factor with levels \code{A} (``theoretical'' departments) 
		or \code{B} (``applied'' departments).}
    \item{\code{yrs.since.phd}}{years since PhD.}
    \item{\code{yrs.service}}{years of service.}
    \item{\code{sex}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{salary}}{nine-month salary, in dollars.}
  }
}

\references{
 Fox J. and Weisberg, S. (2011)
 \emph{An R Companion to Applied Regression}, Second Edition
 Sage.
}

\keyword{datasets}
