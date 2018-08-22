\name{scatterplotMatrix}
\alias{scatterplotMatrix}
\alias{scatterplotMatrix.formula}
\alias{scatterplotMatrix.default}
\alias{spm}
\title{Scatterplot Matrices}
\description{
This function provides a convenient interface to the \code{pairs} function to produce
enhanced scatterplot matrices, including univariate displays on the diagonal and a variety of fitted lines, smoothers, variance functions, and concentration ellipsoids.
\code{spm} is an abbreviation for \code{scatterplotMatrix}.
}
\usage{

scatterplotMatrix(x, ...)

\method{scatterplotMatrix}{formula}(formula, data=NULL, subset, ...)

\method{scatterplotMatrix}{default}(x, smooth = TRUE,
    id = FALSE, legend = TRUE, regLine = TRUE,
    ellipse = FALSE, var.labels = colnames(x), diagonal = TRUE,
    plot.points = TRUE, groups = NULL, by.groups = TRUE,
    use = c("complete.obs", "pairwise.complete.obs"), col =
    carPalette()[-1], pch = 1:n.groups, cex = par("cex"),
    cex.axis = par("cex.axis"), cex.labels = NULL,
    cex.main = par("cex.main"), row1attop = TRUE, ...)

spm(x, ...)
}

\arguments{
  \item{x}{a data matrix or a numeric data frame.}
  \item{formula}{a one-sided \dQuote{model} formula, of the form
    \code{ ~ x1 + x2 + ... + xk} or \code{ ~ x1 + x2 + ... + xk | z} where \code{z}
    evaluates to a factor or other variable to divide the data into groups.}
  \item{data}{for \code{scatterplotMatrix.formula},
    a data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{smooth}{specifies a nonparametric estimate of the mean or median
    function of the vertical axis variable given the
    horizontal axis variable and optionally a nonparametric estimate of the spread or variance function.  If
    \code{smooth=FALSE} neither function is drawn.  If \code{smooth=TRUE}, then both the mean function
    and variance funtions are drawn for ungrouped data, and the mean function only is drawn for grouped
    data.  The default smoother is \code{\link{loessLine}}, which uses the \code{\link{loess}} function from
    the \code{stats} package.  This smoother is fast and reliable.  See the details below for changing
    the smoother, line type, width and color, of the added lines, and adding arguments for the smoother.}
  \item{id}{controls point identification; if \code{FALSE} (the default), no points are identified;
    can be a list of named arguments to the \code{\link{showLabels}} function;
    \code{TRUE} is equivalent to \code{list(method="mahal", n=2, cex=1, location="lr")},
    which identifies the 2 points (in each group, if \code{by.groups=TRUE}) with the largest Mahalanobis distances from the center
    of the data; \code{list(method="identify")} for interactive point identification is not allowed.}
  \item{legend}{controls placement of a legend if the plot is drawn by groups; if \code{FALSE}, the legend
    is suppressed. Can be a list with the named element\code{coords} specifying the position of the legend
    in any form acceptable to the \code{\link{legend}} function;
    \code{TRUE} (the default) is equivalent to \code{list(coords=NULL)}, for which placement will vary
    by the the value of the \code{diagonal} argument---e.g., \code{"topright"} for \code{diagonal=TRUE}.}
   \item{regLine}{controls adding a fitted regression line to each plot, or to each group of points
    if \code{by.groups=TRUE}.  If \code{regLine=FALSE}, no line
    is drawn.  This argument can also be a list with named list, with default \code{regLine=TRUE} equivalent
    to \code{regLine = list(method=lm, lty=1, lwd=2, col=col[1])} specifying the name of the function that
    computes the line, with line type 1 (solid) of relative line width 2 and the color equal to the first
    value in the argument \code{col}. Setting \code{method=MASS::rlm} would fit using a robust regression.}
  \item{ellipse}{controls plotting data-concentration ellipses. If \code{FALSE} (the default), no
    ellipses are plotted.  Can be a list of named values giving \code{levels}, a vector of one or more
    bivariate-normal probability-contour levels at which to
    plot the ellipses;  \code{robust}, a logical value determing whether to use
    the \code{\link{cov.trob}} function in the \pkg{MASS} package
    to calculate the center and covariance matrix for the data ellipses; and \code{fill} and \code{fill.alpha},
    which control whether the ellipse is filled and the transparency of the fill. \code{TRUE} is equivalent
    to \code{list(levels=c(.5, .95), robust=TRUE, fill=TRUE, fill.alpha=0.2)}.}
  \item{var.labels}{variable labels (for the diagonal of the plot).}
  \item{diagonal}{contents of the diagonal panels of the plot. If \code{diagonal=TRUE} adaptive kernel density
    estimates are plotted, separately for each group if grouping is present.  \code{diagonal=FALSE} suppresses
    the diagonal entries.  See details below for other choices for the diagonal.}
  \item{plot.points}{if \code{TRUE} the points are plotted in each
    off-diagonal panel.}
  \item{groups}{a factor or other variable dividing the data into groups; groups are
    plotted with different colors and plotting characters.}
  \item{by.groups}{if \code{TRUE}, the default, regression lines and smooths are fit by groups.}
  \item{use}{if \code{"complete.obs"} (the default), cases with missing data are omitted; if \code{"pairwise.complete.obs"), all valid cases are used
    in each panel of the plot.}}
  \item{pch}{plotting characters for points; default is the plotting characters in
    order (see \code{\link[graphics]{par}}).}
  \item{col}{colors for points; the default is \code{\link{carPalette}} starting at the second color. The color of
    the \code{regLine} and \code{smooth} are the same as for points but can be changed using the the
    \code{regLine} and \code{smooth} arguments.}
  \item{cex}{relative size of plotted points}
  \item{cex.axis}{relative size of axis labels}
  \item{cex.labels}{relative size of labels on the diagonal}
  \item{cex.main}{relative size of the main title, if any}
  \item{row1attop}{If \code{TRUE} (the default) the first row is at the top, as in a matrix, as
  	opposed to at the bottom, as in graph (argument suggested by Richard Heiberger).}
  \item{...}{arguments to pass down.}
}

\details{
   Many arguments to \code{scatterplotMatrix} were changed in version 3 of \pkg{car}, to simplify use of
   this function.

   The \code{smooth} argument is usually either set to \code{TRUE} or \code{FALSE} to draw, or omit,
   the smoother.  Alternatively \code{smooth} can be set to a list of arguments.  The default behavior of
   \code{smooth=TRUE} is equivalent to \code{smooth=list(smoother=loessLine, spread=TRUE, lty.smooth=1, lwd.smooth=1.5, lty.spread=3, lwd.spread=1)}, specifying the smoother to be used, including the spread or variance smooth,
   and the line widths and types for the curves.  You can also specify the colors you want to use for the mean and variance smooths with the arguments \code{col.smooth} and \code{col.spread}. Alternative smoothers are \code{gamline} which uses the
   \code{\link{gam}} function from the \pkg{mgcv} package, and \code{quantregLine} which uses quantile regression to
   estimate the median and quartile functions using \code{\link{rqss}} from the \pkg{quantreg} package.  All of these
   smoothers have one or more arguments described on their help pages, and these arguments can be added to the
   \code{smooth} argument; for example, \code{smooth = list(span=1/2)} would use the default
   \code{loessLine} smoother,
   include the variance smooth, and change the value of the smoothing parameter to 1/2.  For \code{loessLine}
   and \code{gamLine} the variance smooth is estimated by separately
   smoothing the squared positive and negative
   residuals from the mean smooth, using the same type of smoother.  The displayed curves are equal to
   the mean smooth plus the square root of the fit to the positive squared residuals, and the mean fit minus
   the square root of the smooth of the negative squared residuals.  The lines therefore represent the
   comnditional variabiliity at each value on the horizontal axis.  Because smoothing is done separately for
   positive and negative residuals, the variation shown will generally not be symmetric about the fitted mean
   function.  For the \code{quantregLine} method, the center estimates the median for each value on the
   horizontal axis, and the spread estimates the lower and upper quartiles of the estimated conditional
   distribution for each value of the horizontal axis.

      The sub-arguments \code{spread}, \code{lty.spread} and \code{col.spread} of the \code{smooth} argument are equivalent to the newer \code{var}, \code{col.var} and \code{lty.var}, respectively, recognizing that the spread is a measuure of conditional variability.

   By default the diagonal argument is used to draw kernel density estimates of the
   variables by setting \code{diagonal=TRUE}, which is equivalent to setting \code{diagonal =
   list(method="adaptiveDensity", bw=bw.nrd0, adjust=1, kernel=dnorm, na.rm=TRUE)}.  The additional arguments
   shown are descibed at \code{\link{adaptiveKernel}}.  The other methods avaliable, with their default
   arguments, are \code{diagonal=list(method="density", bw="nrd0", adjust=1, kernel="gaussian", na.rm=TRUE)}
   which uses \code{\link{density}} for nonadaptive kernel density estimation; \code{diagonal=list(method
   ="histogram", breaks="FD")}
   which uses \code{\link{hist}} for drawing a histogram that ignores grouping, if present;
   \code{diagonal=list(method="boxplot")} with no additional arguments which draws (parallel) boxplots;
   \code{diagonal=list(method="qqplot")} with no additional arguments which draws a normal QQ plot; and
   \code{diagonal=list(method="oned")} with no additional arguments which draws a rug plot tilted to the
   diagonal, as suggested by Richard Heiberger.

   Earlier versions of \code{scatterplotMatrix} included arguments \code{transform} and \code{family} to estimate power transformations using the \code{\link{powerTransform}} function before drawing the plot.  The same functionality can be achieved by calling \code{powerTransform} directly to estimate a transformation, saving the transformed variables, and then plotting.
}

\value{
  \code{NULL}, returned invisibly. This function is used for its side effect: producing
  a plot.  If point identification is used, a vector of identified points is returned.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\seealso{\code{\link[graphics]{pairs}}, \code{\link{scatterplot}},
  \code{\link{dataEllipse}}, \code{\link{powerTransform}},
  \code{\link{bcPower}}, \code{\link{yjPower}}, \code{\link[MASS]{cov.trob}},
  \code{\link{showLabels}}, \code{\link{ScatterplotSmoothers}}.}

\examples{
scatterplotMatrix(~ income + education + prestige | type, data=Duncan)
scatterplotMatrix(~ income + education + prestige | type, data=Duncan,
    regLine=FALSE, smooth=list(spread=FALSE))
scatterplotMatrix(~ income + education + prestige,
    data=Duncan, id=TRUE, smooth=list(method=gamLine))
}

\keyword{hplot}
