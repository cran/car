\name{scatterplot}
\alias{scatterplot}
\alias{scatterplot.formula}
\alias{scatterplot.default}
\alias{sp}

\title{Enhanced Scatterplots with Marginal Boxplots, Point Marking, Smoothers, and More}
\description{
This function uses basic R graphics to draw a two-dimensional scatterplot, with options to allow for plot enhancements that are often helpful with regression problems.  Enhancements include adding marginal boxplots, estimated mean and variance functions using either parametric or nonparametric methods, point identification, jittering, setting characteristics of points and lines like color, size and symbol, marking points and fitting lines conditional on a grouping variable, and other enhancements.
\code{sp} is an abbreviation for \code{scatterplot}.
}
\usage{
scatterplot(x, ...)

\method{scatterplot}{formula}(formula, data, subset, xlab, ylab, id=FALSE,
    legend=TRUE, ...)

\method{scatterplot}{default}(x, y, boxplots=if (by.groups) "" else "xy",
        regLine=TRUE, legend=TRUE, id=FALSE, ellipse=FALSE, grid=TRUE,
        smooth=TRUE,
        groups, by.groups=!missing(groups),
        xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
        log="", jitter=list(), cex=par("cex"),
        col=carPalette()[-1], pch=1:n.groups,
        reset.par=TRUE, ...)

sp(x, ...)
}

\arguments{


  \item{x}{vector of horizontal coordinates (or first argument of generic function).}
  \item{y}{vector of vertical coordinates.}
  \item{formula}{a model formula, of the form \code{y ~ x} or, if plotting by
    groups, \code{y ~ x | z}, where \code{z} evaluates to a factor
    or other variable dividing the data into groups. If \code{x} is a factor, then parallel boxplots
    are produced using the \code{\link{Boxplot}} function.}
  \item{data}{data frame within which to evaluate the formula.}
  \item{subset}{expression defining a subset of observations.}
  \item{boxplots}{if \code{"x"} a marginal boxplot for the horizontal \code{x}-axis is drawn below the plot;
    if \code{"y"} a marginal boxplot for vertical \code{y}-axis is drawn to the left of the plot;
    if \code{"xy"} both marginal boxplots are drawn; set to \code{""} or \code{FALSE} to
    suppress both boxplots.}
  \item{regLine}{controls adding a fitted regression line to the plot.  if
       \code{regLine=FALSE}, no line is drawn.  If \code{TRUE}, the default, an OLS
       line is fit.  This argument can also be a list.  The default of \code{TRUE} is
       equivalent to \code{regLine=list(method=lm, lty=1, lwd=2, col=col)}, which specifies
       using the \code{lm} function to estimate the fitted line, to draw a solid line
       (\code{lty=1}) of width 2 times the nominal width (\code{lwd=2}) in the color given by
       the first element of the \code{col} argument described below.
       }
  \item{legend}{when the plot is drawn by groups and \code{legend=TRUE}, controls placement
       and properties of a
       legend; if \code{FALSE}, the legend is suppressed. Can be a list of
       named arguments, as follows: \code{title} for the legend; \code{inset}, giving space
      as a proportion of the axes to offset the legend from the axes; \code{coords}
      specifying the position of the legend in any form acceptable to the
      \code{\link{legend}} function or, if not given, the legend is placed \emph{above}
      the plot in the upper margin; \code{columns} for the legend, determined automatically
      to prefer a horizontal layout if not given explicitly; \code{cex} giving the
      relative size of the legend symbols and text. \code{TRUE} (the default) is equivalent to
      \code{list(title=deparse(substitute(groups)), inset=0.02, cex=1)}.}
  \item{id}{controls point identification; if \code{FALSE} (the default), no points are
      identified; can be a list of named arguments to the \code{\link{showLabels}} function;
      \code{TRUE} is equivalent to
      \code{list(method="mahal", n=2, cex=1, col=carPalette()[-1], location="lr")},
      which identifies the 2 points (in each group) with the largest Mahalanobis distances
      from the center of the data.  See \code{\link{showLabels}} for a description of the
      other arguments.  The default behavior of \code{id} is not the same in all graphics
      functions in \pkg{car}, as the \code{method} used depends on the type of plot.}
  \item{ellipse}{controls plotting data-concentration ellipses. If \code{FALSE}
      (the default), no ellipses are plotted.  Can be a list of named values giving
      \code{levels}, a vector of one or more bivariate-normal probability-contour levels at
      which to plot the ellipses; \code{robust}, a logical value determing whether to use
      the \code{\link[MASS]{cov.trob}} function in the \pkg{MASS} package to calculate the center
      and covariance matrix for the data ellipses; and \code{fill} and \code{fill.alpha},
      which control whether the ellipse is filled and the transparency of the fill.
      \code{TRUE} is equivalent to
      \code{list(levels=c(.5, .95), robust=TRUE, fill=TRUE, fill.alpha=0.2)}.}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the graph}
  \item{smooth}{specifies a nonparametric estimate of the mean or median
    function of the vertical axis variable given the
    horizontal axis variable and optionally a nonparametric estimate of the conditional variance.  If
    \code{smooth=FALSE} neither function is drawn.  If \code{smooth=TRUE}, then both the mean function
    and variance funtions are drawn for ungrouped data, and the mean function only is drawn
    for grouped
    data.  The default smoother is \code{\link{loessLine}}, which uses the
    \code{\link{loess}} function from
    the \pkg{stats} package.  This smoother is fast and reliable.  See the details below
    for changing
    the smoother, line type, width and color, of the added lines, and adding arguments
    for the smoother.}
  \item{groups}{a factor or other variable dividing the data into groups; groups are
    plotted with different colors, plotting characters, fits, and smooths.  Using this
    argument is equivalent to specifying the grouping variable in the formula.}
  \item{by.groups}{if \code{TRUE} (the default when there are groups), regression lines are fit by groups.}
  \item{xlab}{label for horizontal axis.}
  \item{ylab}{label for vertical axis.}
  \item{log}{same as the \code{log} argument to \code{\link{plot}}, to produce log axes.}
  \item{jitter}{a list with elements \code{x} or \code{y} or both, specifying jitter factors
    for the horizontal and vertical coordinates of the points in the scatterplot. The
    \code{\link{jitter}} function is used to randomly perturb the points; specifying a
    factor of \code{1} produces the default jitter.
    Fitted lines are unaffected by the jitter.}
  \item{col}{with no grouping, this specifies a color for plotted points;
    with grouping, this argument  should be a vector
    of colors of length at least equal to the number of groups. The default is
    value returned by \code{\link{carPalette}[-1]}.}
  \item{pch}{plotting characters for points; default is the plotting characters in
    order (see \code{\link{par}}).}
  \item{cex}{sets the size of plotting characters, with \code{cex=1} the standard size. You can use \code{cex = 0} to suppress the plotting of points if all you want to show are other graphical features, such as data ellipses, regression lines, smooths, etc.  You can also
    set the sizes of other elements with the arguments \code{cex.axis}, \code{cex.lab}, \code{cex.main},
    and \code{cex.sub}.  See \code{\link{par}}.}
  \item{reset.par}{if \code{TRUE} (the default) then plotting parameters are reset to their previous values
    when \code{scatterplot} exits; if \code{FALSE} then the \code{mar} and \code{mfcol} parameters are
    altered for the current plotting device. Set to \code{FALSE} if you want to add graphical elements
    (such as lines) to the plot.}
  \item{\dots}{other arguments passed down and to \code{plot}.  For example, the argument \code{las} sets
    the style of the axis labels, and \code{xlim} and \code{ylim} set the limits on the horizontal and
    vertical axes, respectively; see \code{\link{par}}.}
}

\details{
   Many arguments to \code{scatterplot} were changed in version 3 of \pkg{car} to simplify use of this function.

The \code{smooth} argument is used to control adding smooth curves to the plot to estimate the conditional center of the vertical axis variable given the horizontal axis variable, and also the conditional variability.  Setting \code{smooth=FALSE} omits all smoothers, while \code{smooth=TRUE}, the default, includes default smoothers.   Alternatively \code{smooth} can be set to a list of subarguments that provide finer control over the smoothing.  

The default behavior of \code{smooth=TRUE} is equivalent to \code{smooth=list(smoother=loessLine, var=!by.groups, lty.var=2, lty.var=4, style="filled", alpha=0.15, border=TRUE, vertical=TRUE)}, specifying the default \code{\link{loessLine}} smoother for the conditional mean smooth and variance smooth.  The color of the smooths is the same of the color of the points, but this can be changed with the arguments \code{col.smooth} and \code{col.var}.
   
Additional available smoothers are \code{\link{gamLine}} which uses the \code{\link[mgcv]{gam}} function and \code{\link{quantregLine}} which uses quantile regression to estimate the median and quartile functions using \code{\link[quantreg]{rqss}}.  All of these smoothers have one or more arguments described on their help pages, and these arguments can be added to the \code{smooth} argument; for example, \code{smooth = list(span=1/2)} would use the default \code{loessLine} smoother, include the variance smooth, and change the value of the smoothing parameter to 1/2.  

For \code{loessLine} and \code{gamLine} the variance smooth is estimated by separately smoothing the squared positive and negative residuals from the mean smooth, using the same type of smoother.  The displayed curves are equal to the mean smooth plus the square root of the fit to the positive squared residuals, and the mean fit minus the square root of the smooth of the negative squared residuals.  The lines therefore represent the comnditional variabiliity at each value on the horizontal axis.  Because smoothing is done separately for positive and negative residuals, the variation shown will generally not be symmetric about the fitted mean function.  For the \code{quantregLine} method, the center estimates the conditional median, and the variability estimates the lower and upper quartiles of the estimated conditional distribution.

The default depection of the variance functions is via a shaded envelope between the upper and lower estimate of variability.  setting the subarguement \code{style="lines"} will display only the boundaries of this region, and \code{style="none"} suppresses variance smoothing.
  
For  \code{style="filled"} several subarguments modify the appearance of the region: \code{alpha} is a number between 0 and 1 that specifies opacity with defualt 0.15, \code{border}, \code{TRUE} or \code{FALSE} specifies a border for the envelope, and \code{vertical} either \code{TRUE} or \code{FALSE}, modifies the behavior of the envelope at the edges of the graph.

The sub-arguments \code{spread}, \code{lty.spread} and \code{col.spread} of the \code{smooth} argument are equivalent to the newer \code{var}, \code{col.var} and \code{lty.var}, respectively, recognizing that the spread is a measuure of conditional variability.
}

\value{
  If points are identified, their labels are returned; otherwise \code{NULL} is returned invisibly.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\references{
  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\seealso{\code{\link{boxplot}},
  \code{\link{jitter}}, \code{\link{legend}},
  \code{\link{scatterplotMatrix}}, \code{\link{dataEllipse}}, \code{\link{Boxplot}},
  \code{\link[MASS]{cov.trob}},
  \code{\link{showLabels}}, \code{\link{ScatterplotSmoothers}}.}

\examples{
scatterplot(prestige ~ income, data=Prestige, ellipse=TRUE, 
  smooth=list(style="lines"))

scatterplot(prestige ~ income, data=Prestige, 
  smooth=list(smoother=quantregLine))
  
scatterplot(prestige ~ income, data=Prestige, 
  smooth=list(smoother=quantregLine, border="FALSE"))

# use quantile regression for median and quartile fits
scatterplot(prestige ~ income | type, data=Prestige,
  smooth=list(smoother=quantregLine, var=TRUE, span=1, lwd=4, lwd.var=2))

scatterplot(prestige ~ income | type, data=Prestige, 
  legend=list(coords="topleft"))

scatterplot(vocabulary ~ education, jitter=list(x=1, y=1),
            data=Vocab, smooth=FALSE, lwd=3)

scatterplot(infantMortality ~ ppgdp, log="xy", data=UN, id=list(n=5))

scatterplot(income ~ type, data=Prestige)

\dontrun{ # interactive point identification
    # remember to exit from point-identification mode
    scatterplot(infantMortality ~ ppgdp, id=list(method="identify"), data=UN)
}
}


\keyword{hplot}
