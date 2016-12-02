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

\method{Boxplot}{default}(y, g, labels, id.method = c("y", "identify", "none"), 
    id.n=10, xlab, ylab, ...)

\method{Boxplot}{formula}(formula, data = NULL, subset, na.action = NULL, labels., 
    id.method = c("y", "identify", "none"), xlab, ylab, ...)

\method{Boxplot}{list}(y, xlab="", ylab="", ...)

\method{Boxplot}{data.frame}(y, labels=rownames(y), ...)

\method{Boxplot}{matrix}(y, ...)
}

\arguments{
  \item{y}{a numeric variable for which the boxplot is to be constructed; a list of numeric
    variables, each element of which will be treated as a group; a numeric data frame or a
    numeric matrix, each of whose columns will be treated as a group.}
  \item{g}{a grouping variable, usually a factor, for constructing parallel boxplots.}
  \item{labels, labels.}{point labels; if not specified, \code{Boxplot} will use the row names of the \code{data}
    argument, if one is given, or observation numbers, or row names if \code{y} is a data frame or matrix (that has row names).}
  \item{id.method}{if \code{"y"} (the default), all outlying points are labeled; if \code{"identify"}, points
    may be labeled interactive; if \code{"none"}, no point identification is performed.}
  \item{id.n}{up to \code{id.n} high outliers and low outliers will be identified in each group, (default, 10).}
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
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
}

\seealso{
\code{\link{boxplot}}
}

\examples{
Boxplot(~income, data=Prestige, id.n=Inf) # identify all outliers
Boxplot(income ~ type, data=Prestige)
Boxplot(income ~ type, data=Prestige, at=c(1, 3, 2))
Boxplot(k5 + k618 ~ lfp*wc, data=Mroz)
with(Prestige, Boxplot(income, labels=rownames(Prestige)))
with(Prestige, Boxplot(income, type, labels=rownames(Prestige)))
Boxplot(scale(Prestige[, 1:4]))
}

\keyword{hplot}
