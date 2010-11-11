\name{symbox}
\alias{symbox}
\alias{symbox.formula}
\alias{symbox.default}

\title{Boxplots for transformations to symmetry}

\description{\code{symbox} first transforms \code{x} to each of a series
    of selected powers, with each transformation standardized to mean 0 and
    standard deviation 1. The results are then displayed side-by-side in boxplots, permiting
    a visual assessment of which power makes the distribution reasonably symmetric.}

\usage{
symbox(x, ...)
\method{symbox}{formula}(formula, data=NULL, subset, na.action=NULL, ylab,  ...)
\method{symbox}{default}(x, powers = c(-1, -0.5, 0, 0.5, 1), start=0, 
	trans=bcPower, xlab="Powers", ylab, ...)
}

\arguments{
    \item{x}{a numeric vector.}
    \item{formula}{a one-sided formula specifying a single numeric variable.}
    \item{data, subset, na.action}{as for statistical modeling functions (see, e.g., \code{\link{lm}}).}
    \item{xlab, ylab}{axis labels; if \code{ylab} is missing, a label will be supplied.}
    \item{powers}{a vector of selected powers to which \code{x} is to be raised. 
        For meaningful comparison of powers,  \code{1} should be included
        in the vector of powers.}
    \item{start}{a constant to be added to \code{x}.}
    \item{trans}{a transformation function whose first argument is a numeric vector and whose second
    	argument is a transformation parameter, given by the \code{powers} argument; the default
    	is \code{\link{bcPower}}, and another possibility is \code{\link{yjPower}}.}
    \item{\ldots}{arguments to be passed down.}
    }

\value{as returned by \code{boxplot}.}

\author{Gregor Gorjanc, John Fox \email{jfox@mcmaster.ca}, and Sanford Weisberg.}

\references{

	Fox, J. and Weisberg, S. (2011)
    \emph{An R Companion to Applied Regression}, Second Edition. Sage.
	
}
    
\seealso{\code{\link[graphics]{boxplot}}, \code{\link[MASS]{boxcox}}, \code{\link{bcPower}}, \code{\link{yjPower}}}

\examples{
symbox(~ income, data=Prestige)
}

\keyword{hplot}
