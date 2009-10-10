\name{recode}
\alias{recode}
\title{Recode a Variable}
\description{
 Recodes a numeric vector, character vector, or factor 
 according to simple recode specifications.
}
\usage{
recode(var, recodes, as.factor.result, levels)
}

\arguments{
  \item{var}{numeric vector, character vector, or factor.}
  \item{recodes}{character string of recode specifications:
    see below.}
  \item{as.factor.result}{return a factor; default is \code{TRUE} if 
    \code{var} is a factor, \code{FALSE} otherwise.}
  \item{levels}{an optional argument specifying the order of the
    levels in the returned factor; the default is to use the sort order
    of the level names.}
}
\details{
  Recode specifications appear in a character string, separated by
  semicolons (see the examples below), of the form \code{input=output}.
  If an input value satisfies more than one specification, 
  then the first (from left to right) applies.
  If no specification is satisfied, then the input value is carried
  over to the result. \code{NA} is allowed on input and output.
  Several recode specifications are supported:
  \describe{
    \item{single value}{For example, \code{0=NA}.}
    \item{vector of values}{For example, \code{c(7,8,9)='high'}.}
    \item{range of values}{For example, \code{7:9='C'}. The special values \code{lo}
      and \code{hi} may appear in a range. For example, \code{lo:10=1}.}
    \item{\code{else}}{everything that does not fit a previous specification.
      For example, \code{else=NA}. Note that \code{else} matches \emph{all} otherwise
      unspecified values on input, including \code{NA}.}
    }
  If all of the output values are numeric, and if \code{as.factor.result} is 
  \code{FALSE}, then a numeric result is returned.
}

\value{
  a recoded vector of the same length as \code{var}; if \code{var} is a factor,
  then so is the result.
}
 
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[base]{cut}}, \code{\link[base]{factor}}}

\examples{
x<-rep(1:3,3)
x
## [1] 1 2 3 1 2 3 1 2 3
recode(x, "c(1,2)='A'; else='B'")
## [1] "A" "A" "B" "A" "A" "B" "A" "A" "B"
recode(x, "1:2='A'; 3='B'")
## [1] "A" "A" "B" "A" "A" "B" "A" "A" "B"
}

\keyword{manip}
