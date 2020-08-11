\name{recode}
\alias{recode}
\alias{Recode}
\title{Recode a Variable}
\description{
 Recodes a numeric vector, character vector, or factor 
 according to simple recode specifications. \code{Recode}
 is an alias for \code{recode} that avoids name clashes
 with packages, such as \pkg{Hmisc}, that have a \code{recode} function.
}
\usage{
recode(var, recodes, as.factor, as.numeric=TRUE, levels)

Recode(...)
}

\arguments{
  \item{var}{numeric vector, character vector, or factor.}
  \item{recodes}{character string of recode specifications:
    see below.}
  \item{as.factor}{return a factor; default is \code{TRUE} if 
    \code{var} is a factor, \code{FALSE} otherwise.}
  \item{as.numeric}{if \code{TRUE} (the default), and \code{as.factor} is \code{FALSE}, 
  	then the result will be coerced to numeric
  	if all values in the result are numerals---i.e., represent numbers.}
  \item{levels}{an optional argument specifying the order of the
    levels in the returned factor; the default is to use the sort order
    of the level names.}
  \item{...}{arguments to be passed to \code{recode}.}
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
      and \code{hi} may appear in a range. For example, \code{lo:10=1}. \emph{Note:} \code{:} is
      \emph{not} the R sequence operator.  In addition, you may not use \code{:} with the \code{c} function within a recode specification, so for example \code{c(1, 3, 5:7)} will cause an error.}
    \item{\code{else}}{everything that does not fit a previous specification.
      For example, \code{else=NA}. Note that \code{else} matches \emph{all} otherwise
      unspecified values on input, including \code{NA}.}
    }
  Character data and factor levels on the left-hand side of a recode specification must be quoted. Thus,
  e.g., \code{c(a, b, c) = 'low'} is not allowed, and should be \code{c('a', 'b', 'c') = 'low'}. 
  Similarly, the colon is reserved for numeric data, and, e.g., \code{c('a':'c') = 'low'} is not allowed. 
  If the \code{var} argument is a character variable with (some) values that are composed of numerals, or a factor
  with (some) levels that are numerals (e.g., \code{'12'} or \code{'-2'}), then these too must be quoted 
  and cannot be used with colons (e.g., \code{'15':'19' = '15 to 19'} is not allowed, and could be
  specified as \code{c('15', '16', '17', '18', '19') = '15 to 19'}, assuming that all values are 
  the character representation of whole numbers).
  
  If all of the output values are numeric, and if \code{as.factor} is 
  \code{FALSE}, then a numeric result is returned; if \code{var} is a factor,
  then by default so is the result.
}

\value{
  a recoded vector of the same length as \code{var}.
}
 
\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\seealso{\code{\link{cut}}, \code{\link{factor}}}

\examples{
x<-rep(1:3,3)
x
## [1] 1 2 3 1 2 3 1 2 3
recode(x, "c(1,2)='A'; 
	else='B'")
## [1] "A" "A" "B" "A" "A" "B" "A" "A" "B"
Recode(x, "1:2='A'; 3='B'")
## [1] "A" "A" "B" "A" "A" "B" "A" "A" "B"
}

\keyword{manip}
