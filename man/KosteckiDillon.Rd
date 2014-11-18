\name{KosteckiDillon}
\alias{KosteckiDillon}
\docType{data}
\title{Treatment of Migraine Headaches}

\description{
Subset of data on migraine treatments collected by Tammy Kostecki-Dillon. 
}

\usage{KosteckiDillon}

\format{
  A data frame with 4152 observations on 133 subjects for the following 9 variables.
  \describe{
    \item{\code{id}}{Patient id.}
    \item{\code{time}}{time in days relative to the onset of treatment, which occurs at time 0.}
    \item{\code{dos}}{time in days from the start of the study, January 1 of the first year of the study.}
    \item{\code{hatype}}{a factor with levels \code{Aura} \code{Mixed} \code{No Aura}, the type of
        migraine experienced by a subject.}
    \item{\code{age}}{at onset of treatment, in years.}
    \item{\code{airq}}{a measure of air quality.}
    \item{\code{medication}}{a factor with levels \code{none} \code{reduced} \code{continuing},
        representing subjects who discontinued their medication, who continued but at a reduced
        dose, or who continued at the previous dose.}
    \item{\code{headache}}{a factor with levels \code{no} \code{yes}.}
    \item{\code{sex}}{a factor with levels \code{female} \code{male}.}
  }
}

\details{
The data consist of headache logs kept by 133 patients in a treatment program in which bio-feedback 
was used to attempt to reduce migraine frequency and severity.
Patients entered the program at different times over a period of about 3 years. 
Patients were encouraged to begin their logs four weeks before the onset of treatment 
and to continue for one month afterwards, 
but only 55 patients have data preceding the onset of treatment.
}

\source{
Personal communication from Georges Monette (and adapted from his description of the data).
}

\references{
Kostecki-Dillon, T.,  Monette, G., and Wong, P. (1999). 
Pine trees, comas, and migraines.
\emph{York University Institute for Social Research Newsletter},
14:2.
}

\examples{
summary(KosteckiDillon)
}
\keyword{datasets}
