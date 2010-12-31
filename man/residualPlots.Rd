\name{residualPlots}
\alias{residualPlots}
\alias{residualPlots.default}
\alias{residualPlots.lm}
\alias{residualPlots.glm}
\alias{residualPlot}
\alias{residualPlot.default}
\alias{residualPlot.lm}
\alias{residualPlot.glm}
\alias{residCurvTest}
\alias{residCurvTest.lm}
\alias{residCurvTest.glm}
\alias{tukeyNonaddTest}
\title{Residual Plots and Curvature Tests for Linear Model Fits}
\description{
Plots the residuals versus each term in a mean function and versus
fitted values.  Also computes a curvature test for each of the plots
by adding a quadratic term and testing the quadratic to be zero.  This
is Tukey's test for nonadditivity when plotting against fitted values.
}
\usage{
### This is a generic function with only one required argument:

residualPlots (model, ...)

\S3method{residualPlots}{default}(model, terms = ~., layout = NULL, ask, 
                 main = "", fitted = TRUE, AsIs=FALSE, plot = TRUE, 
                 tests = TRUE, ...)

\S3method{residualPlots}{lm}(model, ...)

\S3method{residualPlots}{glm}(model, ...) 

### residualPlots calls residualPlot, so these arguments can be
### used with either function

residualPlot(model, ...)  

\S3method{residualPlot}{default}(model, variable = "fitted", type = "pearson", 
                 plot = TRUE,     
                 quadratic = TRUE, 
                 smooth = FALSE, span = 1/2, smooth.lwd=lwd, smooth.lty=lty,
                 smooth.col=col.lines, 
                 labels,
                 id.method = "y",
                 id.n = if(id.method[1]=="identify") Inf else 0,
                 id.cex=1, id.col=palette()[1],
                 col = palette()[1], col.lines = palette()[2], 
                 xlab, ylab, lwd = 1, lty=1, grid=TRUE, ...)   
                 
\S3method{residualPlot}{lm}(model, ...) 
  
\S3method{residualPlot}{glm}(model, variable = "fitted", type = "pearson", 
                 plot = TRUE, quadratic = FALSE, smooth = TRUE, ...)  
}
\arguments{
  \item{model}{ 
A regression object.
}
  \item{terms}{
A one-sided formula that specifies a subset of the predictors. One
residual plot is drawn for each specified.  The default
\code{~.} is to plot against all predictors.  For example, the
specification \code{terms = ~.-X3} would plot against all predictors
except for \code{X3}.  Some nonstandard predictors, such as B-splines, are 
skipped.
}
  \item{layout}{ 
If set to a value like \code{c(1, 1)} or \code{c(4, 3)}, the layout
of the graph will have this many rows and columns.  If not set, the program
will select an appropriate layout.  If the number of graphs exceed nine, you 
must select the layout yourself, or you will get a maximum of nine per page.
}
  \item{ask}{
If \code{TRUE}, ask the user before drawing the next plot; if \code{FALSE}, don't
ask.
}
  \item{main}{
Main title for the graphs.  The default is \code{main=""} for no title.
}
  \item{fitted}{
If \code{TRUE}, the default, include the plot against fitted values.
}
  \item{AsIs}{
If \code{FALSE}, the default, terms that use the \dQuote{as-is} function \code{I}
are skipped; if \code{TRUE}, they are included.
}
  \item{plot}{
If \code{TRUE}, draw the plot(s).
}
  \item{tests}{
If \code{TRUE}, display the curvature tests.
}
  \item{...}{
Additional arguments passed to \code{residualPlot} and then to
\code{plot}.
}
  \item{variable}{
Quoted variable name for the horizontal axis, or 
\code{"fitted"} to plot versus fitted values.
}
  \item{type}{
Type of residuals to be used.  Pearson residuals are 
appropriate for \code{lm} objects since these are equivalent to ordinary residuals
with ols and correctly weighted residuals with wls.  Any quoted string that
is an appropriate value of the \code{type} argument to
\code{\link{residuals.lm}} or \code{"rstudent"} or \code{"rstandard"} for 
Studentized or standardized residuals.
}
  \item{quadratic}{
if \code{TRUE}, fits the quadratic regression of the
vertical axis on the horizontal axis and displays a lack of fit test.  Default
is \code{TRUE} for \code{lm} and \code{}FALSE for \code{glm}.  
}
  \item{smooth}{
if TRUE fits a loess smooth using the settings given below.  Defaults to FALSE
for lm objects and TRUE for glm objects.
}
  \item{span, smooth.lwd, smooth.lty, smooth.col}{
Should a lowess smooth be added to the figure?  The \code{span} is the smoothing
parameter for lowess, \code{smooth.lwd}, \code{smooth.lty}, and \code{smooth.col} 
are, respectively, the width, type, and color of the line drawn on the plot. 
}
  \item{id.method,labels,id.n,id.cex,id.col}{Arguments for the labelling of
    points.  The default is \code{id.n=0} for labeling no points.  See
    \code{\link{showLabels}} for details of these arguments.
}
  \item{col}{
default color for points
}
  \item{col.lines}{
default color for lines
}
  \item{xlab}{
X-axis label.  If not specified, a useful label is constructed by
the function.
}
  \item{ylab}{
Y-axis label.  If not specified, a useful label is constructed by
the function.
}
  \item{lwd}{
line width for lines.
}
  \item{lty}{
line type for quadratic.
}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
  }
\details{
\code{residualPlots} draws one or more residuals plots depending on the
value of the \code{terms} and \code{fitted} arguments.  If \code{terms = ~ .},
the default, then a plot is produced of residuals versus each first-order 
term in the formula used to create the model.  Interaction terms, spline terms,
and polynomial terms of more than one predictor are 
skipped.  In addition terms that use the \dQuote{as-is} function, e.g., \code{I(X^2)},
will also be skipped unless you set the argument \code{AsIs=TRUE}.  A plot of 
residuals versus fitted values is also included unless \code{fitted=FALSE}. 

In addition to plots, a table of curvature tests is displayed. For plots
against a term in the model formula, say \code{X1}, the test displayed is
the t-test for for \code{I(X^2)} in the fit of \code{update, model, ~. + I(X^2))}.
Econometricians call this a specification test.  For factors, the displayed
plot is a boxplot, and no curvature test is computed.
For fitted values, the test is Tukey's one-degree-of-freedom test for 
nonadditivity.  You can suppress the tests with the argument \code{tests=FALSE}.
  
\code{residualPlot}, which is called by \code{residualPlots},
should be viewed as an internal function, and is included here to display its
arguments, which can be used with \code{residualPlots} as well.  The
\code{residualPlot} function returns the curvature test as an invisible result.
  
\code{residCurvTest} computes the curvature test only.  For any factors a
boxplot will be drawn.  For any polynomials, plots are against the linear term.
Other non-standard predictors like B-splines are skipped.
}
\value{
For \code{lm} objects,
returns a data.frame with one row for each plot drawn, one column for 
the curvature test statistic, and a second column for the corresponding
p-value.  This function is used primarily for its side effect of drawing
residual plots. 
}
\references{

Fox, J. and Weisberg, S. (2011)
\emph{An R Companion to Applied Regression}, Second Edition.  Sage.

Weisberg, S. (2005) \emph{Applied
Linear Regression}, Third Edition, Wiley, Chapter 8}
\author{Sanford Weisberg, \email{sandy@stat.umn.edu}}

\seealso{See Also \code{\link{lm}}, \code{\link{identify}},
\code{\link{showLabels}} }


\examples{
residualPlots(lm(longley))
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
