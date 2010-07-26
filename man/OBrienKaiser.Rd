\name{OBrienKaiser}
\alias{OBrienKaiser}
\docType{data}
\title{O'Brien and Kaiser's Repeated-Measures Data}
\description{
  These contrived repeated-measures data are taken from 
  O'Brien and Kaiser (1985). The data are from an imaginary study in which
  16 female and male subjects, who are divided into three treatments, are measured
  at a pretest, postest, and a follow-up session; during each session, they are
  measured at five occasions at intervals of one hour. The design, therefore, has
  two between-subject and two within-subject factors.
  
  The contrasts for the \code{treatment} factor are set to \eqn{-2, 1, 1} and
  \eqn{0, -1, 1}. The contrasts for the \code{gender} factor are set to
  \code{contr.sum}.
}

\usage{OBrienKaiser}

\format{
  A data frame with 16 observations on the following 17 variables.
  \describe{
    \item{\code{treatment}}{a factor with levels \code{control} \code{A} \code{B}}
    \item{\code{gender}}{a factor with levels \code{F} \code{M}}
    \item{\code{pre.1}}{pretest, hour 1}
    \item{\code{pre.2}}{pretest, hour 2}
    \item{\code{pre.3}}{pretest, hour 3}
    \item{\code{pre.4}}{pretest, hour 4}
    \item{\code{pre.5}}{pretest, hour 5}
    \item{\code{post.1}}{posttest, hour 1}
    \item{\code{post.2}}{posttest, hour 2}
    \item{\code{post.3}}{posttest, hour 3}
    \item{\code{post.4}}{posttest, hour 4}
    \item{\code{post.5}}{posttest, hour 5}
    \item{\code{fup.1}}{follow-up, hour 1}
    \item{\code{fup.2}}{follow-up, hour 2}
    \item{\code{fup.3}}{follow-up, hour 3}
    \item{\code{fup.4}}{follow-up, hour 4}
    \item{\code{fup.5}}{follow-up, hour 5}
  }
}

\source{
  O'Brien, R. G., and Kaiser, M. K. (1985)
  MANOVA method for analyzing repeated measures designs: An extensive primer.
  \emph{Psychological Bulletin} \bold{97}, 316--333, Table 7.
}

\examples{
OBrienKaiser
contrasts(OBrienKaiser$treatment)
contrasts(OBrienKaiser$gender)
}
\keyword{datasets}
