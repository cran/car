\name{WeightLoss}
\Rdversion{1.1}
\alias{WeightLoss}
\docType{data}
\title{
Weight Loss Data
}
\description{
Contrived data on weight loss and self esteem over three months,
for three groups of individuals: Control, Diet and Diet + Exercise.
The data constitute a double-multivariate design.
}
\usage{WeightLoss}
\format{
  A data frame with 34 observations on the following 7 variables.
  \describe{
    \item{\code{group}}{a factor with levels \code{Control} \code{Diet} \code{DietEx}.}
    \item{\code{wl1}}{Weight loss at 1 month}
    \item{\code{wl2}}{Weight loss at 2 months}
    \item{\code{wl3}}{Weight loss at 3 months}
    \item{\code{se1}}{Self esteem at 1 month}
    \item{\code{se2}}{Self esteem at 2 months}
    \item{\code{se3}}{Self esteem at 3 months}
  }
}
\details{
Helmert contrasts are assigned to \code{group}, comparing \code{Control} vs. (\code{Diet} \code{DietEx})
and \code{Diet} vs. \code{DietEx}.
}
\source{
Originally taken from \url{http://www.csun.edu/~ata20315/psy524/main.htm},
but modified slightly. Courtesy of Michael Friendly.
}

\keyword{datasets}
