\name{spreadLevelPlot}
\alias{spreadLevelPlot}
\alias{slp}
\alias{spreadLevelPlot.formula}
\alias{spreadLevelPlot.default}
\alias{spreadLevelPlot.lm}
\alias{print.spreadLevelPlot}

\title{Spread-Level Plots}
\description{
  Creates plots for examining the possible dependence of
  spread on level, or an extension of these plots to the studentized residuals
  from linear models.
}
\usage{
spreadLevelPlot(x, ...)

slp(...)

\method{spreadLevelPlot}{formula}(x, data=NULL, subset, na.action, 
    main=paste("Spread-Level Plot for", varnames[response], 
    "by", varnames[-response]), ...)

\method{spreadLevelPlot}{default}(x, by, robust.line=TRUE, 
	start=0, xlab="Median", ylab="Hinge-Spread", 
	point.labels=TRUE, las=par("las"),
	main=paste("Spread-Level Plot for", deparse(substitute(x)), 
		"by", deparse(substitute(by))), 
	col=carPalette()[1], col.lines=carPalette()[2],
    pch=1, lwd=2, grid=TRUE, ...)

\method{spreadLevelPlot}{lm}(x, robust.line=TRUE, 
	xlab="Fitted Values", ylab="Absolute Studentized Residuals", las=par("las"),
	main=paste("Spread-Level Plot for\n", deparse(substitute(x))),
	pch=1, col=carPalette()[1], col.lines=carPalette()[2:3], lwd=2, grid=TRUE, 
    id=FALSE, smooth=TRUE, ...)
	  
\method{print}{spreadLevelPlot}(x, ...)
}

\arguments{
  \item{x}{a formula of the form \code{y ~ x}, where \code{y} is a numeric vector
    and \code{x} is a factor, or an \code{lm} object to be plotted; alternatively a numeric vector.}
  \item{data}{an optional data frame containing the variables to be plotted. 
    By default the variables are taken from the environment from which 
    \code{spreadLevelPlot} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  \item{na.action}{a function that indicates what should happen when the data contain \code{NA}s. 
    The default is set by the \code{na.action} setting of \code{options}.}
  \item{by}{a factor, numeric vector, or character vector defining groups.}
  \item{robust.line}{if \code{TRUE} a robust line is fit using the \code{rlm} function in
    the \code{MASS} package; if \code{FALSE} a line is fit using \code{lm}.}
  \item{start}{add the constant \code{start} to each data value.}
  \item{main}{title for the plot.}
  \item{xlab}{label for horizontal axis.}
  \item{ylab}{label for vertical axis.}
  \item{point.labels}{if \code{TRUE} label the points in the plot with group names.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{col}{color for points; the default is the first entry
    in the current \pkg{car} palette (see \code{\link{carPalette}}
    and \code{\link[graphics]{par}}).}
  \item{col.lines}{for the default method, the line color, defaulting to the second entry in the \pkg{car} color palette; 
    for the \code{"lm"} method, a vector of two colors for, respectively, the fitted straight line and a nonparametric regression smooth line,
    default to the second and third entries in the \pkg{car} color palette.}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
    \item{id}{controls point identification; if \code{FALSE} (the default), no points are identified;
        can be a list of named arguments to the \code{\link{showLabels}} function;
        \code{TRUE} is equivalent to \code{list(method=list("x", "y"), n=2, 
        cex=1, col=carPalette()[1], location="lr")},
        which identifies the 2 points the most extreme horizontal ("X", absolute studentized residual) values and the 2 points
        with the most extreme horizontal ("Y", fitted values) values.}
    \item{smooth}{specifies the smoother to be used along with its arguments; if \code{FALSE}, no smoother is shown;
        can be a list giving the smoother function and its named arguments; \code{TRUE}, the default, is equivalent to
        \code{list(smoother=loessLine)}. See \code{\link{ScatterplotSmoothers}} for the smoothers supplied by the
        \pkg{car} package and their arguments.}
  \item{\dots}{arguments passed to plotting functions.}
}
\details{
  Except for linear models, computes the statistics for, and plots, a Tukey spread-level plot
  of log(hinge-spread) vs. log(median) for the groups; fits a line to the plot; and calculates a
  spread-stabilizing transformation from the slope of the line. 
  
  For linear models, plots log(abs(studentized residuals) vs. log(fitted values).  Point labeling was added in November, 2016.
  
  The function \code{slp} is an abbreviation for \code{spreadLevelPlot}.
}
\value{
  An object of class \code{spreadLevelPlot} containing:
  \item{Statistics}{a matrix with the lower-hinge, median, upper-hinge, and hinge-spread
    for each group. (Not for an \code{lm} object.)}
  \item{PowerTransformation}{spread-stabilizing power transformation, calculated as \eqn{1 - slope}
    of the line fit to the plot.}
}
\references{ 
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
  
  Hoaglin, D. C., Mosteller, F. and Tukey, J. W. (Eds.) (1983)
  \emph{Understanding Robust and Exploratory Data Analysis.} Wiley.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{hccm}}, \code{\link{ncvTest}} }

\examples{
spreadLevelPlot(interlocks + 1 ~ nation, data=Ornstein)
slp(lm(interlocks + 1 ~ assets + sector + nation, data=Ornstein))
}

\keyword{hplot}
\keyword{regression}
