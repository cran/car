\name{strings2factors}
\alias{strings2factors}
\alias{strings2factors.data.frame}

\title{
Convert Character-String Variables in a Data Frame to Factors
}

\description{
Converts the character variables (or a subset of these variables) in a data frame to factors, with optional control of the order of the resulting factor levels.
}

\usage{
strings2factors(object, which, not, exclude.unique, levels, verbose, ...)
\method{strings2factors}{data.frame}(object, which, not, 
    exclude.unique=TRUE, levels=list(), verbose=TRUE, ...)
}

\arguments{
  \item{object}{a data frame or an object inheriting from the \code{"data.frame"} class.}
  \item{which}{an optional character vector of names or column numbers of the character variables to be converted to factors; if absent, \emph{all} character variables will be converted, except as excluded by the \code{not} and \code{exclude.unique} arguments (see below).}
  \item{not}{an optional character vector of names or column numbers of character variables \emph{not} to be converted to factors.}
    \item{exclude.unique}{if \code{TRUE} (the default), character variables all of whose values are unique (i.e., all different from each other) are not converted to factors. Such variables, which would have as many levels as there are cases, are typically case identifiers and not categorical variables. If \code{FALSE}, character variables all of whose values are unique are converted to factors with a warning.}
  \item{levels}{an optional named list, each element of which is a character vector of levels of the corresponding factor. This argument allows you to control the order of levels of the factor; if omitted, or if a particular factor is omitted from the list, the levels will be in the default alphabetic order.}
  \item{verbose}{if \code{TRUE} (the default), the names of the character variables that were converted to factors are printed on the console.}
  \item{\dots}{not used.}
}

\value{a data frame with (some) character variables replaced by corresponding factors.}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{factor}}, \code{\link{data.frame}}}

\examples{
M <- Moore # from the carData package
M$partner <- as.character(Moore$partner.status)
M$fcat <- as.character(Moore$fcategory)
M$names <- rownames(M) # values are unique
str(M)
str(strings2factors(M))
str(strings2factors(M,
  levels=list(partner=c("low", "high"), fcat=c("low", "medium", "high"))))
str(strings2factors(M, which="partner", levels=list(partner=c("low", "high"))))
str(strings2factors(M, not="partner", exclude.unique=FALSE))
}

\keyword{manip}
