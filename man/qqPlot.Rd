\name{qqPlot}
\alias{qqPlot}
\alias{qqp}
\alias{qqPlot.default}
\alias{qqPlot.formula}
\alias{qqPlot.lm}
\alias{qqPlot.glm}

\title{Quantile-Comparison Plot}
\description{
 Plots empirical quantiles of a variable, or of studentized residuals from
 a linear model, against theoretical quantiles of a comparison distribution.  Includes
 options not avaiable in the \code{\link[stats]{qqnorm}} function.
}
\usage{
qqPlot(x, ...)

qqp(...)

\method{qqPlot}{default}(x, distribution="norm", groups, layout,
    ylim=range(x, na.rm=TRUE), ylab=deparse(substitute(x)),
    xlab=paste(distribution, "quantiles"), glab=deparse(substitute(groups)),
    main=NULL, las=par("las"),
    envelope=.95, col=carPalette()[1], col.lines=carPalette()[2],
    lwd=2, pch=1, cex=par("cex"),
    line=c("quartiles", "robust", "none"), id=TRUE, grid=TRUE, ...)

\method{qqPlot}{formula}(formula, data, subset, id=TRUE, ylab, glab, ...)

\method{qqPlot}{lm}(x, xlab=paste(distribution, "Quantiles"),
    ylab=paste("Studentized Residuals(",
                deparse(substitute(x)), ")", sep=""),
    main=NULL, distribution=c("t", "norm"),
    line=c("robust", "quartiles", "none"), las=par("las"),
    simulate=TRUE, envelope=.95,  reps=100,
    col=carPalette()[1], col.lines=carPalette()[2], lwd=2, pch=1, cex=par("cex"),
    id=TRUE, grid=TRUE, ...)
}

\arguments{
  \item{x}{vector of numeric values or \code{lm} object.}
    \item{distribution}{root name of comparison distribution -- e.g., \code{"norm"} for the
    normal distribution; \code{t} for the t-distribution.}
  \item{groups}{an optional factor; if specified, a QQ plot will be drawn for \code{x}
    within each level of \code{groups}.}
  \item{layout}{a 2-vector with the number of rows and columns for plotting by
    groups -- for example \code{c(1, 3)} for 1 row and 3 columns; if omitted, the
    number of rows and columns will be selected automatically; the specified number
    of rows and columns must be sufficient to accomodate the number of groups; ignored
    if there is no grouping factor.}
  \item{formula}{one-sided formula specifying a single variable to be plotted or a two-sided formula of
    the form \code{variable ~ factor}, where a QQ plot will be drawn for \code{variable} within each
    level of \code{factor}.}
  \item{data}{optional data frame within which to evaluage the formula.}
  \item{subset}{optional subset expression to select cases to plot.}
  \item{ylim}{limits for vertical axis; defaults to the range of \code{x}. If plotting by groups, a common
    y-axis is used for all groups.}
  \item{ylab}{label for vertical (empirical quantiles) axis.}
  \item{xlab}{label for horizontal (comparison quantiles) axis.}
  \item{glab}{label for the grouping variable.}
  \item{main}{label for plot.}
  \item{envelope}{confidence level for point-wise confidence envelope, or
    \code{FALSE} for no envelope.}
  \item{las}{if \code{0}, ticks labels are drawn parallel to the
    axis; set to \code{1} for horizontal labels (see \code{\link[graphics]{par}}).}
  \item{col}{color for points; the default is the \emph{first} entry
    in the current \pkg{car} palette (see \code{\link{carPalette}}
    and \code{\link[graphics]{par}}).}
  \item{col.lines}{color for lines; the default is the \emph{second} entry
    in the current \pkg{car} palette.}
  \item{pch}{plotting character for points; default is \code{1}
    (a circle, see \code{\link[graphics]{par}}).}
  \item{cex}{factor for expanding the size of plotted symbols; the default is
    \code{1}.}
  \item{id}{controls point identification; if \code{FALSE}, no points are identified;
    can be a list of named arguments to the \code{\link{showLabels}} function;
    \code{TRUE} is equivalent to \code{list(method="y", n=2, cex=1, col=carPalette()[1], location="lr")},
    which identifies the 2 points with the 2 points with the most extreme
    verical values --- studentized residuals for the \code{"lm"} method. Points labels are by default
    taken from the names of the variable being plotted is any, else case indices are used.  Unlike most graphical functions in \pkg{car}, the default is \code{id=TRUE} to include point identification.}
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
  of an independent random sample from the comparison distribution (see Fox, 2016).
  Studentized residuals from linear models are plotted against the appropriate t-distribution with a point-wise
  confidence envelope computed by default by a parametric bootstrap,
  as described by Atkinson (1985).
  The function \code{qqp} is an abbreviation for \code{qqPlot}.
}
\value{
  These functions return the labels of identified points, unless a grouping factor is employed,
  in which case \code{NULL} is returned invisibly.
}
\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.

  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.

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

qqPlot(~ income, data=Prestige, subset = type == "prof")
qqPlot(income ~ type, data=Prestige, layout=c(1, 3))

qqPlot(lm(prestige ~ income + education + type, data=Duncan),
	envelope=.99)
}

\keyword{distribution}
\keyword{univar}
\keyword{regression}
