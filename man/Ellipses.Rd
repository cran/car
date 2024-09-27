\name{Ellipses}
\alias{ellipse}
\alias{dataEllipse}
\alias{dataEllipse.formula}
\alias{dataEllipse.default}
\alias{confidenceEllipse}
\alias{confidenceEllipse.default}
\alias{confidenceEllipse.lm}
\alias{confidenceEllipse.glm}
\alias{confidenceEllipse.mlm}
\alias{confidenceEllipses}
\alias{confidenceEllipses.default}
\alias{confidenceEllipses.mlm}

\title{Ellipses, Data Ellipses, and Confidence Ellipses}
\description{
  These functions draw ellipses, including data ellipses, and
  confidence ellipses for linear, generalized linear, and possibly other models.
}
\usage{
ellipse(center, shape, radius, log="", center.pch=19, center.cex=1.5, 
  segments=51, draw=TRUE, add=draw, xlab="", ylab="", 
  col=carPalette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, grid=TRUE, ...)

dataEllipse(x, ...)

\method{dataEllipse}{formula}(formula, data, subset, weights, xlab, ylab, id=FALSE, ...)

\method{dataEllipse}{default}(x, y, groups, 
  group.labels=group.levels, ellipse.label, 
  weights, log="", levels=c(0.5, 0.95), center.pch=19, 
  center.cex=1.5, draw=TRUE, plot.points=draw, add=!plot.points, segments=51, 
  robust=FALSE, 
  xlab=deparse(substitute(x)), 
  ylab=deparse(substitute(y)), 
  col=if (missing(groups)) carPalette()[1:2] else carPalette()[1:length(group.levels)],
  pch=if (missing(groups)) 1 else seq(group.levels),
  lwd=2, fill=FALSE, fill.alpha=0.3, grid=TRUE, id=FALSE, 
  label.pos=NULL, label.cex=1.25, label.xpd=FALSE, ...)

confidenceEllipse(model, ...)

\method{confidenceEllipse}{lm}(model, which.coef, vcov.=vcov, 
  L, levels=0.95, Scheffe=FALSE, dfn,
  center.pch=19, center.cex=1.5, segments=51, xlab, ylab, 
  col=carPalette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, draw=TRUE, add=!draw,
  grid=TRUE, ...)
  
\method{confidenceEllipse}{glm}(model, chisq, ...)

\method{confidenceEllipse}{mlm}(model, xlab, ylab, which.coef=1:2, ...)

\method{confidenceEllipse}{default}(model, which.coef, vcov.=vcov,
  L, levels=0.95, Scheffe=FALSE,  dfn,
  center.pch=19, center.cex=1.5, segments=51, xlab, ylab, 
  col=carPalette()[2], lwd=2, fill=FALSE, fill.alpha=0.3, draw=TRUE, add=!draw, 
  grid=TRUE, ...)
  
confidenceEllipses(model, ...)

\method{confidenceEllipses}{default}(model, coefnames,  main, grid=TRUE, ...)

\method{confidenceEllipses}{mlm}(model, coefnames, main, ...)
}

\arguments{
  \item{center}{2-element vector with coordinates of center of ellipse.}
  \item{shape}{\eqn{2\times 2}{2 * 2} shape (or covariance) matrix.}
  \item{radius}{radius of circle generating the ellipse.}
  \item{log}{when an ellipse is to be added to an existing plot, indicates
    whether computations were on logged values and to be plotted on logged
    axes; \code{"x"} if the x-axis is logged, \code{"y"} if the y-axis is
    logged, and \code{"xy"} or \code{"yx"} if both axes are logged. The
    default is \code{""}, indicating that neither axis is logged.}
  \item{center.pch}{character for plotting ellipse center; if \code{FALSE} or \code{NULL} the center point isn't plotted.}
  \item{center.cex}{relative size of character for plotting ellipse center.}
  \item{segments}{number of line-segments used to draw ellipse.}
  \item{draw}{if \code{TRUE} produce graphical output; if \code{FALSE}, only invisibly return coordinates of ellipse(s).}
  \item{add}{if \code{TRUE} add ellipse(s) to current plot.}
  \item{xlab}{label for horizontal axis.}
  \item{ylab}{label for vertical axis.}
  \item{x}{a numeric vector, or (if \code{y} is missing) a 2-column numeric matrix.}
  \item{y}{a numeric vector, of the same length as \code{x}.}
  \item{groups}{optional: a factor to divide the data into groups; a separate ellipse will be plotted for each group (level of the factor).}
  \item{formula}{for the formula method, a model formula, of the form \code{y ~ x} or, if plotting by
    groups, \code{y ~ x | z}, where \code{z} evaluates to a factor
    or other variable dividing the data into groups.}
  \item{data}{for the formula method, a data frame within which to evaluate the formula.}
  \item{subset}{for the formula method, expression defining a subset of observations.}
  \item{group.labels}{labels to be plotted for the groups; by default, the levels of the \code{groups} factor.}
  \item{ellipse.label}{a label for the ellipse(s) or a vector of labels; if several ellipses are drawn and just one label is given,
    then that label will be repeated. The default is not to label the ellipses.}
  \item{label.pos}{position(s) of labels for the ellipses, a vector of integers (in 0:4) or character strings (in c("center", "bottom", "left", "top", "right"), or in c("C", "S", "W", "N", "E") used in labeling ellipses, recycled as necessary. Values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the max/min coordinates of the ellipse; the value 0 specifies the centroid of the ellipse object. The default, \code{label.pos=NULL} uses the correlation of the ellipse to determine "top" (r>=0) or "bottom" (r<0)}.
  \item{label.cex}{character size for ellipse labels.}
  \item{label.xpd}{allow labels to be drawn outside of plotting region (default \code{FALSE}).}
  \item{weights}{a numeric vector of weights, of the same length as \code{x} and \code{y} to be used  by \code{\link{cov.wt}}
    or \code{\link[MASS]{cov.trob}} in computing a weighted covariance matrix; if absent, weights of \code{1} are used.}
  \item{plot.points}{if \code{FALSE} data ellipses are drawn,
    but points are not plotted.}
  \item{levels}{draw elliptical contours at these (normal) probability or confidence levels.}
  \item{robust}{if \code{TRUE} use the \code{cov.trob} function in the \pkg{MASS} package
    to calculate the center and covariance matrix for the data ellipse.}
  \item{model}{a model object produced by \code{lm} or \code{glm}.}
  \item{which.coef}{2-element vector giving indices of coefficients to plot; 
    if missing, the first two coefficients (disregarding the regression constant)
    will be selected.}
  \item{vcov.}{a coefficient-covariance matrix or a function (such as \code{\link{hccm}}) to 
    compute the coefficent-covariance
    matrix from \code{model}; the default is the \code{vcov} function.
    
    Note that arguments supplied to \code{...} are \emph{not} passed to \code{vcov.} when it's a function; in this case either use an anonymous function in which the 
    additional arguments are set, or supply the coefficient covariance matrix directly.
    }
  \item{L}{As an alternative to selecting coefficients to plot, a transformation matrix can be specified to compute two
    linear combinations of the coefficients; if the \code{L} matrix is given, it takes precedence over the \code{which.coef}
    argument. \code{L} should have two rows and as many columns as there are coefficients. It can be given directly as a
    numeric matrix, or specified by a pair of character-valued expressions, in the same manner as for the 
    \code{link{linearHypothesis}} function, but with no right-hand side.}
  \item{Scheffe}{if \code{TRUE} scale the ellipse so that its projections onto the
    axes give Scheffe confidence intervals for the coefficients.}
  \item{dfn}{``numerator'' degrees of freedom (or just degrees of freedom for a GLM) for
  	drawing the confidence ellipse. Defaults to the number of coefficients in the model (disregarding the constant) if 
  	\code{Scheffe} is \code{TRUE}, or to \code{2} otherwise; selecting \code{dfn = 1} will 
  	draw the ``confidence-interval generating'' ellipse, with projections on the axes
  	corresponding to individual confidence intervals with the stated level of coverage.}
  \item{chisq}{if \code{TRUE}, the confidence ellipse for the coefficients in a generalized linear model are
    based on the chisquare statistic, if \code{FALSE} on the $F$-statistic. This corresponds to using the default
    and linear-model methods respectively.}
  \item{col}{color for lines and ellipse center; the default is the \emph{second} entry
    in the current \pkg{car} palette (see \code{\link{carPalette}}
    and \code{\link{par}}). For \code{dataEllipse}, two colors can be given, in
    which case the first is for plotted points and the second for lines and the ellipse center;
    if ellipses are plotted for \code{groups}, then this is a vector of colors for the groups.}
  \item{pch}{for \code{dataEllipse} this is the plotting character (default, symbol \code{1}, a hollow circle)
    to use for the points; if ellipses are plotted by \code{groups}, then this a vector of plotting characters,
    with consecutive symbols starting with \code{1} as the default.}
  \item{lwd}{line width; default is \code{2} (see \code{\link{par}}).}
  \item{fill}{fill the ellipse with translucent color \code{col} (default, \code{FALSE})?}
  \item{fill.alpha}{transparency of fill (default = \code{0.3}).}
  \item{\dots}{other plotting parameters to be passed to \code{plot} and
    \code{line}.}
  \item{id}{controls point identification; if \code{FALSE} (the default), no points are identified;
    can be a list of named arguments to the \code{\link{showLabels}} function;
    \code{TRUE} is equivalent to \code{list(method="mahal", n=2, cex=1, col=carPalette()[1], location="lr")}
    (with the default \code{col} actually dependent on the number of groups),
    which identifies the 2 points with the largest Mahalanobis distances from the center of the data.}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
  \item{coefnames}{character vector of coefficient names to use to label the diagonal of the pairwise confidence ellipse matrix plotted by \code{confidenceEllipses}; defaults to the names of the coefficients in the model.}
  \item{main}{title for matrix of pairwise confidence ellipses.}
}
\details{
  The ellipse is computed by suitably transforming a unit circle.
  
  \code{dataEllipse} superimposes the normal-probability contours over a scatterplot
  of the data.
  
  \code{confidenceEllipses} plots a matrix of all pairwise confidence ellipses; each panel of the matrix is created by \code{confidenceEllipse}.
}
\value{
   These functions are mainly used for their side effect of producing plots.  For
   greater flexibility (e.g., adding plot annotations), however, \code{ellipse} returns invisibly the (x, y) coordinates of the calculated ellipse.
   \code{dataEllipse} and \code{confidenceEllipse} return invisibly the  coordinates of one or more ellipses, in the latter instance a list named by
   \code{levels}; \code{confidenceEllipses} invisibly returns \code{NULL}.
}

\references{ 
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
  
  Monette, G. (1990)
  Geometry of multiple regression and 3D graphics.
  In Fox, J. and Long, J. S. (Eds.) 
  \emph{Modern Methods of Data Analysis.} Sage.
 }
 
\author{Georges Monette, John Fox \email{jfox@mcmaster.ca}, and Michael Friendly.}

\seealso{\code{\link[MASS]{cov.trob}}, \code{\link{cov.wt}}, \code{\link{linearHypothesis}}.}
  
\examples{
dataEllipse(Duncan$income, Duncan$education, levels=0.1*1:9, 
    ellipse.label=0.1*1:9, lty=2, fill=TRUE, fill.alpha=0.1)

# formula method
dataEllipse(education ~ income | type, data=Duncan,
            levels = 0.68)

dataEllipse(education ~ income | type, data=Duncan, 
            levels=0.68, label.xpd=TRUE)

# vary label positions
dataEllipse(education ~ income | type, data=Duncan,
            levels = 0.68,
            fill = TRUE, fill.alpha = 0.1,
            label.pos = c(3, 2, 4),
            xlim = c(0, 100),
            ylim = c(2, 120))

    
confidenceEllipse(lm(prestige ~ income + education, data=Duncan), Scheffe=TRUE)

confidenceEllipse(lm(prestige ~ income + education, data=Duncan), vcov.=hccm)

confidenceEllipse(lm(prestige ~ income + education, data=Duncan), 
	L=c("income + education", "income - education"))
	
confidenceEllipses(lm(prestige ~ income + education + type, data=Duncan),
  fill=TRUE)
cov2cor(vcov(lm(prestige ~ income + education + type, 
  data=Duncan))) # correlations among coefficients

wts <- rep(1, nrow(Duncan))
wts[c(6, 16)] <- 0 # delete Minister, Conductor
with(Duncan, {
	dataEllipse(income, prestige, levels=0.68)
	dataEllipse(income, prestige, levels=0.68, robust=TRUE, 
	    plot.points=FALSE, col="green3")
	dataEllipse(income, prestige, weights=wts, levels=0.68, 
	    plot.points=FALSE, col="brown")
	dataEllipse(income, prestige, weights=wts, robust=TRUE, levels=0.68, 
		plot.points=FALSE, col="blue")
	})
    
with(Prestige, dataEllipse(income, education, type, 
    id=list(n=2, labels=rownames(Prestige)), pch=15:17,
    xlim=c(0, 25000), center.pch="+",
    group.labels=c("Blue Collar", "Professional", "White Collar"),
    ylim=c(5, 20), level=.95, fill=TRUE, fill.alpha=0.1))
}

\keyword{hplot}
\keyword{aplot}
