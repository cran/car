\name{Boxplot}
\alias{Boxplot}
\alias{Boxplot.default}
\alias{Boxplot.formula}
\alias{Boxplot.list}
\alias{Boxplot.data.frame}
\alias{Boxplot.matrix}

\title{
Boxplots With Point Identification
}

\description{
\code{Boxplot} is a wrapper for the standard \R{} \code{\link{boxplot}} function, providing point identification, 
axis labels, and a formula interface for boxplots without a grouping variable.
}

\usage{
Boxplot(y, ...)

\method{Boxplot}{default}(y, g, id=TRUE, xlab, ylab, ...)

\method{Boxplot}{formula}(formula, data=NULL, subset, na.action=NULL, 
    id=TRUE, xlab, ylab, ...)

\method{Boxplot}{list}(y, xlab="", ylab="", ...)

\method{Boxplot}{data.frame}(y, id=TRUE, ...)

\method{Boxplot}{matrix}(y, ...)
}

\arguments{
  \item{y}{a numeric variable for which the boxplot is to be constructed; a list of numeric
    variables, each element of which will be treated as a group; a numeric data frame or a
    numeric matrix, each of whose columns will be treated as a group.}
  \item{g}{a grouping variable, usually a factor, for constructing parallel boxplots.}
  \item{id}{a list of named elements giving one or more specifications for labels of individual points ("outliers"):
    \code{n}, the maximum number of points to label (default 10); \code{location}, \code{"lr"} (left or right) of points or 
    \code{"avoid"} to try to avoid overplotting; \code{method}, one of \code{"y"} (automatic, the default), \code{"identify"} (interactive), or \code{"none"}; \code{col} for labels (default is the first color in \code{carPalette()} ); and \code{cex} size of labels (default is \code{1}).
    Can be \code{FALSE} to suppress point identification or \code{TRUE} (the default) to use all defaults. This is similar to how \code{\link{showLabels}} handles point labels
    for other functions in the \pkg{car} package, except that the usual default is \code{id=FALSE}.}
  \item{xlab, ylab}{text labels for the horizontal and vertical axes; if missing, \code{Boxplot} will use the
    variable names, or, in the case of a list, data frame, or matrix, empty labels.}
  \item{formula}{a `model' formula, of the form \code{~ y} to produce a boxplot for the variable \code{y}, or
    of the form \code{y ~ g}, \code{y ~ g1*g2*...}, or \code{y ~ g1 + g2 + ...} to 
    produce parallel boxplots for \code{y} within levels of the grouping variable(s)
    \code{g}, etc., usually factors.}
  \item{data, subset, na.action}{as for statistical modeling functions (see, e.g., \code{\link{lm}}).}
  \item{\dots}{further arguments, such as \code{at}, to be passed to \code{\link{boxplot}}.}
}

\author{John Fox \email{jfox@mcmaster.ca}, with a contribution from Steve Ellison 
    to handle \code{at} argument (see \code{\link{boxplot}}).}

\references{
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\seealso{
\code{\link{boxplot}}
}

\examples{
Boxplot(~income, data=Prestige, id=list(n=Inf)) # identify all outliers
Boxplot(income ~ type, data=Prestige)
Boxplot(income ~ type, data=Prestige, at=c(1, 3, 2))
Boxplot(k5 + k618 ~ lfp*wc, data=Mroz)
with(Prestige, Boxplot(income, id=list(labels=rownames(Prestige))))
with(Prestige, Boxplot(income, type, id=list(labels=rownames(Prestige))))
Boxplot(scale(Prestige[, 1:4]))
}

\keyword{hplot}
