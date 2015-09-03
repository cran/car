\name{Wong}
\alias{Wong}
\docType{data}
\title{
Post-Coma Recovery of IQ
}

\description{
The \code{Wong} data frame has 331 row and 7 columns. The observations are longitudinal
data on recovery of IQ after comas of varying duration for 200 subjects.
}

\usage{Wong}

\format{
  This data frame contains the following columns:
  \describe{
    \item{\code{id}}{patient ID number.}
    \item{\code{days}}{number of days post coma at which IQs were measured.}
    \item{\code{duration}}{duration of the coma in days.}
    \item{\code{sex}}{a factor with levels \code{Female} and \code{Male}.}
    \item{\code{age}}{in years at the time of injury.}
    \item{\code{piq}}{performance (i.e., mathematical) IQ.}
    \item{\code{viq}}{verbal IQ.}
  }
}

\details{
The data are from Wong, Monette, and Weiner (2001) and are for 200 patients who sustained traumatic brain
injuries resulting in comas of varying duration. After
awakening from their comas, patients were periodically administered a standard IQ test, but the average
number of measurements per patient is small (331/200 = 1.7).
}

\source{
Wong, P. P., Monette, G., and Weiner, N. I. (2001) Mathematical models of cognitive recovery. \emph{Brain Injury}, \bold{15}, 519--530.
}

\references{
Fox, J. (2016) \emph{Applied Regression Analysis and Generalized Linear Models}, Third Edition. Sage.
}

\examples{
summary(Wong)
}

\keyword{datasets}
