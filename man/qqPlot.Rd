\name{qqPlot}
\alias{qqPlot}
\alias{qqp}
\alias{qqPlot.default}
\alias{qqPlot.lm}
\alias{qqPlot.glm}

\title{Quantile-Comparison Plots}
\description{
 Plots empirical quantiles of a variable, or of studentized residuals from
 a linear model, against theoretical quantiles of a comparison distribution.
}
\usage{
qqPlot(x, ...)

qqp(...)

\method{qqPlot}{default}(x, distribution="norm", ...,
	ylab=deparse(substitute(x)), xlab=paste(distribution, "quantiles"), 
	main=NULL, las=par("las"),
	envelope=.95,  
	col=palette()[1], col.lines=palette()[2], lwd=2, pch=1, cex=par("cex"), 
	line=c("quartiles", "robust", "none"),  
	labels = if(!is.null(names(x))) names(x) else seq(along=x),
	id.method = "y", 
	id.n =if(id.method[1]=="identify") Inf else 0,
	id.cex=1, id.col=palette()[1], grid=TRUE)

\method{qqPlot}{lm}(x, xlab=paste(distribution, "Quantiles"),
	ylab=paste("Studentized Residuals(", deparse(substitute(x)), ")",
		sep=""), main=NULL,
	distribution=c("t", "norm"), line=c("robust", "quartiles", "none"),
	las=par("las"), simulate=TRUE, envelope=.95,  
	reps=100, col=palette()[1], col.lines=palette()[2], lwd=2, 
	pch=1, cex=par("cex"),
	labels, id.method = "y", 
	id.n = if(id.method[1]=="identify") Inf else 0,
	id.cex=1, id.col=palette()[1], grid=TRUE, ...)
}

\arguments{
  \item{x}{vector of numeric values or \code{lm} object.}
  \item{distribution}{root name of comparison distribution -- e.g., \code{"norm"} for the
    normal distribution; \code{t} for the t-distribution.}
  \item{ylab}{label for vertical (empirical quantiles) axis.}
  \item{xlab}{label for horizontal (comparison quantiles) axis.}
  \item{main}{label for plot.}
  \item{envelope}{confidence level for point-wise confidence envelope, or 
    \code{FALSE} for no envelope.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{col}{color for points; the default is the \emph{first} entry
    in the current color palette (see \code{\link[grDevices]{palette}}
    and \code{\link[graphics]{par}}).}
  \item{col.lines}{color for lines; the default is the \emph{second} entry
    in the current color palette.}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).}
  \item{cex}{factor for expanding the size of plotted symbols; the default is
    \code{1}.}
  \item{labels}{vector of text strings to be used to identify points, defaults to
    \code{names(x)} or observation numbers if \code{names(x)} is \code{NULL}.}
  \item{id.method}{point identification method.  The default
\code{id.method="y"}
    will identify the \code{id.n} points with the largest value of 
    \code{abs(y-mean(y))}.  See \code{\link{showLabels}} for other options.}
  \item{id.n}{number of points labeled.  If \code{id.n=0}, the default, no
    point identification.} 
  \item{id.cex}{set size of the text for point labels; the default is \code{cex} (which is typically \code{1}).}
  \item{id.col}{color for the point labels.}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{line}{\code{"quartiles"} to pass a line through the quartile-pairs, or
    \code{"robust"} for a robust-regression line; the latter uses the \code{rlm}
    function in the \code{MASS} package. Specifying \code{line = "none"} suppresses the line.}
  \item{simulate}{if \code{TRUE} calculate confidence envelope by parametric bootstrap;
    for \code{lm} object only. The method is due to Atkinson (1985).}
  \item{reps}{integer; number of bootstrap replications for confidence envelope.}
  \item{\dots}{arguments such as \code{df} to be passed to the appropriate quantile function.}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
}

\details{
  Draws theoretical quantile-comparison plots for variables and for studentized residuals
  from a linear model. A comparison line is drawn on the plot either through the quartiles
  of the two distributions, or by robust regression. 
  
  Any distribution for which quantile and
  density functions exist in R (with prefixes \code{q} and \code{d}, respectively) may be used. 
  When plotting a vector, the confidence envelope is based on the SEs of the order statistics 
  of an independent random sample from the comparison distribution (see Fox, 2008).
  Studentized residuals from linear models are plotted against the appropriate t-distribution with a point-wise 
  confidence envelope computed by default by a parametric bootstrap,
  as described by Atkinson (1985).
  The function \code{qqp} is an abbreviation for \code{qqPlot}.
}
\value{
  These functions return the labels of identified points.
}
\references{
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Atkinson, A. C. (1985)
  \emph{Plots, Transformations, and Regression.} Oxford.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[stats]{qqplot}}, \code{\link[stats]{qqnorm}},
 \code{\link[stats]{qqline}}, \code{\link{showLabels}}}

\examples{
x<-rchisq(100, df=2)
qqPlot(x)
qqPlot(x, dist="chisq", df=2)

qqPlot(lm(prestige ~ income + education + type, data=Duncan), 
	envelope=.99)
}

\keyword{distribution}
\keyword{univar}
\keyword{regression}
