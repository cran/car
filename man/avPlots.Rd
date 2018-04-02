\name{avPlots}
\alias{avPlots}
\alias{avp}
\alias{avPlot}
\alias{avPlot.lm}
\alias{avPlot.glm}

\title{Added-Variable Plots}
\description{
  These functions construct added-variable, also called 
  partial-regression, plots for  linear and generalized linear
  models.
}

\usage{
avPlots(model, terms=~., intercept=FALSE, layout=NULL, ask, main, ...)

avp(...)

avPlot(model, ...)

\method{avPlot}{lm}(model, variable,
	id=TRUE, col = carPalette()[1], col.lines = carPalette()[2],
	xlab, ylab, pch = 1, lwd = 2, 
	main=paste("Added-Variable Plot:", variable),
	grid=TRUE,
	ellipse=FALSE,
  marginal.scale=FALSE, ...)

\method{avPlot}{glm}(model, variable, 
	id=TRUE,
	col = carPalette()[1], col.lines = carPalette()[2],
	xlab, ylab, pch = 1, lwd = 2,  type=c("Wang", "Weisberg"), 
	main=paste("Added-Variable Plot:", variable), grid=TRUE,
	ellipse=FALSE, ...)
}
\arguments{
  \item{model}{model object produced by \code{lm} or \code{glm}.}
  \item{terms}{ A one-sided formula that specifies a subset of the predictors. 
    One added-variable plot is drawn for each term.   For example, the
    specification \code{terms = ~.-X3} would plot against all terms
    except for \code{X3}. If this argument is a quoted name of one of the terms, the 
    added-variable plot is drawn for that term only.}
  \item{intercept}{Include the intercept in the plots; default is \code{FALSE}.}
  \item{variable}{A quoted string giving the name of a regressor in the model matrix for the horizontal
    axis.}
  \item{layout}{If set to a value like \code{c(1, 1)} or \code{c(4, 3)}, the layout
    of the graph will have this many rows and columns.  If not set, the program
    will select an appropriate layout.  If the number of graphs exceed nine, you 
    must select the layout yourself, or you will get a maximum of nine per page.
    If \code{layout=NA}, the function does not set the layout and the user can
    use the \code{par} function to control the layout, for example to have 
    plots from two models in the same graphics window.}
  \item{main}{The title of the plot; if missing, one will be supplied.}
  \item{ask}{If \code{TRUE}, ask the user before drawing the next plot; if \code{FALSE} don't ask.}
  \item{\dots}{\code{avPlots} passes these arguments to \code{avPlot}.  
    \code{avPlot} passes them to \code{plot}. }
  \item{id}{controls point identification; if \code{FALSE}, no points are identified;
    can be a list of named arguments to the \code{\link{showLabels}} function;
    \code{TRUE}, the default, is equivalent to 
    \code{list(method=list(abs(residuals(model, type="pearson")), "x"), n=2, cex=1, col=carPalette()[1], location="lr")},
    which identifies the 2 points with the largest residuals and the 2 points with 
    the most extreme horizontal values (i.e., largest partial leverage).}
  \item{col}{color for points; the default is the \emph{second} entry
    in the current \pkg{car} palette (see \code{\link{carPalette}}
    and \code{\link[graphics]{par}}).}
  \item{col.lines}{color for the fitted line.}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{xlab}{x-axis label.  If omitted a label will be constructed.}
  \item{ylab}{y-axis label.  If omitted a label will be constructed.}
  \item{type}{if \code{"Wang"} use the method of Wang (1985);
    if \code{"Weisberg"} use the method in the Arc software associated with
    Cook and Weisberg (1999).}
  \item{grid}{If \code{TRUE}, the default, a light-gray background grid is put on the graph.}
  \item{ellipse}{controls plotting data-concentration ellipses. If \code{FALSE} (the default), no ellipses are plotted.
    Can be a list of named values giving \code{levels}, a vector of one or more bivariate-normal probability-contour levels at which to
    plot the ellipses; and \code{robust}, a logical value determing whether to use the \code{\link{cov.trob}} function in the \pkg{MASS} package
    to calculate the center and covariance matrix for the data ellipses. \code{TRUE} is equivalent to \code{list(levels=c(.5, .95), robust=TRUE)}.}  
  \item{marginal.scale}{Consider an added-variable plot of Y versus X given Z.  If this argument is \code{FALSE} then the limits on the horizontal axis are determined by the range of the residuals from the regression of X on Z and the limits on the vertical axis are determined by the range of the residuals from the regressnio of Y on Z.  If the argument is \code{TRUE}, then the limits on the horizontal axis are determined by the range of X minus it mean, and on the vertical axis by the range of Y minus its means; adjustment is made if necessary to include outliers.  This scaling allows visualization of the correlations between Y and Z and between X and Z.  For example, if the X and Z are highly correlated, then the points will be concentrated on the middle of the plot.} 
}

\details{
  The function intended for direct use is \code{avPlots} (for which \code{avp}
  is an abbreviation). 
}

\value{
  These functions are used for their side effect id producing
  plots, but also invisibly return the coordinates of the plotted points.
}

\references{
  Cook, R. D. and Weisberg, S. (1999)
  \emph{Applied Regression, Including Computing and Graphics.} Wiley.

  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
  
  Wang, P C. (1985)
  Adding a variable in generalized linear models.
  \emph{Technometrics} \bold{27}, 273--276.
  
  Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}}

\seealso{\code{\link{residualPlots}}, \code{\link{crPlots}}, \code{\link{ceresPlots}}, \code{link{dataEllipse}}, \code{\link{showLabels}}, \code{\link{dataEllipse}}.
}

\examples{
avPlots(lm(prestige ~ income + education + type, data=Duncan))

avPlots(glm(partic != "not.work" ~ hincome + children, 
  data=Womenlf, family=binomial), id=FALSE)
  
m1 <- lm(partic ~ tfr + menwage + womwage + debt + parttime, Bfox)
par(mfrow=c(1,3))
# marginal plot, ignoring other predictors:
with(Bfox, dataEllipse(womwage, partic, levels=0.5)) 
abline(lm(partic ~ womwage, Bfox), col="red", lwd=2)
# AV plot, adjusting for others:
avPlots(m1, ~ womwage, ellipse=list(levels=0.5)) 
# AV plot, adjusting and scaling as in marginal plot
avPlots(m1, ~ womwage, marginal.scale=TRUE, ellipse=list(levels=0.5)) 
}

\keyword{hplot}
\keyword{regression}
