\name{car-internal}
\alias{inv}
\alias{has.intercept}
\alias{has.intercept.default}
\alias{term.names}
\alias{term.names.default}
\alias{predictor.names}
\alias{predictor.names.default}
\alias{response.name}
\alias{response.name.default}
\alias{response}
\alias{response.default}
\alias{is.aliased}
\alias{df.terms}
\alias{df.terms.default}
\alias{mfrow}
\alias{nice}
\alias{relatives}

\title{Internal car functions}

\description{
Internal functions for package car.
}

\usage{
df.terms(model, term, ...)
df.terms.default(model, term, ...)
has.intercept(model, ...)
has.intercept.default(model)
inv(x)
is.aliased(model)
mfrow(n, max.plots=0)
nice(x, direction=c("round", "down", "up"))
predictor.names(model, ...)
predictor.names.default(model)
relatives(term, names, factors)
response.name (model, ...) 
response.name.default(model)
response(model, ...)
response.default(model)
term.names(model, ...)
term.names.default(model)
}

\arguments{
  \item{model}{}
  \item{term}{}
  \item{...}{}
  \item{x}{}
  \item{n}{}
  \item{max.plots}{}
  \item{direction}{}
  \item{names}{}
  \item{factors}{}
}

\details{
These functions are not intended to be called by the user.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\keyword{misc}
