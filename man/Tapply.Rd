\name{Tapply}
\alias{Tapply}
\title{Apply a Function to a Variable Within Factor Levels}
\description{
Applies a function, typically to compute a single statistic, like a mean, median, or standard deviation,
within levels of a factor or within combinations of levels of two or more factors to produce a table of statistics. This function provides
a formula interface to the standard R \code{\link{tapply}} function.
}
 
\usage{
Tapply(formula, fun, data, na.action = na.pass, ..., targs = list())
}
\arguments{
  \item{formula}{ a two-sided formula of the form \code{variable ~ factor.1 + factor.2 + ... + factor.n} or,
  equivalently, \code{variable ~ factor.1*factor.2* ... *factor.n}.
  The variables on the right-hand side
  of the formula are normally factors or are otherwise coerced to factors.}
  \item{fun}{a function, like \code{mean}, \code{median}, or \code{sd}, that takes a vector first argument and typically
  returns a single number as its value.}
  \item{data}{an optional data frame within which to find the variable and factor(s).}
  \item{na.action}{a function to handle missing values, as in statistical modeling functions like \code{\link{lm}};
  the default is \code{\link{na.pass}}.}
  \item{\dots}{arguments to pass to the function given in the \code{fun} argument, such as (commonly) 
  \code{na.rm=TRUE}.}
  \item{targs}{a list of optional arguments to pass to \code{\link{tapply}}.}
}
\details{
The function given by \code{fun} is applied to the values of the left-hand-side variable in \code{formula} within
(combination of) levels of the factor(s) given in the right-hand side of \code{formula}, producing a table of statistics.
}

\value{
The object returned by \code{\link{tapply}}, typically simply printed.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
	Fox, J. and Weisberg, S. (2019)
    \emph{An R Companion to Applied Regression}, Third Edition. Sage.
}


\seealso{\code{\link{tapply}}.}

\examples{
Tapply(conformity ~ partner.status + fcategory, mean, data=Moore)
Tapply(conformity ~ partner.status + fcategory, mean, data=Moore, 
    trim=0.2)

Moore[1, 2] <- NA
Tapply(conformity ~ partner.status + fcategory, mean, data=Moore)
Tapply(conformity ~ partner.status + fcategory, mean, data=Moore, 
  na.rm=TRUE)
Tapply(conformity ~ partner.status + fcategory, mean, data=Moore, 
  na.action=na.omit)  # equivalent
remove("Moore")
}
\keyword{misc}
\keyword{manip}
