\name{qq.plot}
\alias{qq.plot}
\alias{qqp}
\alias{qq.plot.default}
\alias{qq.plot.lm}
\alias{qq.plot.glm}

\title{Quantile-Comparison Plots}
\description{
 Plots empirical quantiles of a variable, or of studentized residuals from
 a linear model, against theoretical quantiles of a comparison distribution.
}
\usage{
qq.plot(x, ...)

qqp(...)

\method{qq.plot}{default}(x, distribution="norm", ylab=deparse(substitute(x)),
        xlab=paste(distribution, "quantiles"), main=NULL, las=par("las"),
        envelope=.95, labels=FALSE, col=palette()[2], lwd=2, pch=1, cex=1,
        line=c("quartiles", "robust", "none"), ...)

\method{qq.plot}{lm}(x, main=NULL, xlab=paste(distribution, "Quantiles"), 
  ylab=paste("Studentized Residuals(", deparse(substitute(x)), ")",
    sep = ""), 
  distribution=c("t", "norm"), line=c("quartiles", "robust", "none"), 
  las=par('las'), simulate=FALSE, envelope=0.95, labels=names(rstudent), 
  reps=100, col=palette()[2], lwd=2, pch=1, cex=1, ...)
}
\arguments{
  \item{x}{vector of numeric values or \code{lm} object.}
  \item{distribution}{root name of comparison distribution -- e.g., \code{norm} for the
    normal distribution; \code{t} for the t-distribution.}
  \item{ylab}{label for vertical (empirical quantiles) axis.}
  \item{xlab}{label for horizontal (comparison quantiles) axis.}
  \item{main}{label for plot.}
  \item{envelope}{confidence level for point-wise confidence envelope, or 
    \code{FALSE} for no envelope.}
  \item{labels}{vector of point labels for interactive point identification,
    or \code{FALSE} for no labels.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{col}{color for points and lines; the default is the \emph{second} entry
    in the current color palette (see \code{\link[grDevices]{palette}}
    and \code{\link[graphics]{par}}).}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).}
  \item{cex}{factor for expanding the size of plotted symbols; the default is
    \code{1}.}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).
    Confidence envelopes are drawn at half this line width.}
  \item{line}{\code{"quartiles"} to pass a line through the quartile-pairs, or
    \code{"robust"} for a robust-regression line; the latter uses the \code{rlm}
    function in the \code{MASS} package. Specifying \code{line = "none"} suppresses the line.}
  \item{simulate}{if \code{TRUE} calculate confidence envelope by parametric bootstrap;
    for \code{lm} object only. The method is due to Atkinson (1985).}
  \item{reps}{integer; number of bootstrap replications for confidence envelope.}
  \item{\dots}{arguments such as \code{df} to be passed to the appropriate quantile function.}
}

\details{
  Draws theoretical quantile-comparison plots for variables and for studentized residuals
  from a linear model. A comparison line is drawn on the plot either through the quartiles
  of the two distributions, or by robust regression. 
  
  Any distribution for which quantile and
  density functions exist in R (with prefixes \code{q} and \code{d}, respectively) may be used. 
  Studentized residuals are plotted against the
  appropriate t-distribution.
  
  The function \code{qqp} is an abbreviation for \code{qq.plot}.
}
\value{
  \code{NULL}. These functions are used only for their side effect (to make a graph).
}
\references{
  Fox, J. (1997)
  \emph{Applied Regression, Linear Models, and Related Methods.} Sage.

  Atkinson, A. C. (1985)
  \emph{Plots, Transformations, and Regression.} Oxford.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[stats]{qqplot}}, \code{\link[stats]{qqnorm}},
 \code{\link[stats]{qqline}}}

\examples{
x<-rchisq(100, df=2)
qq.plot(x)
qq.plot(x, dist="chisq", df=2)

qq.plot(lm(interlocks~assets+sector+nation, data=Ornstein), sim=TRUE)
}

\keyword{distribution}
\keyword{univar}
\keyword{regression}
