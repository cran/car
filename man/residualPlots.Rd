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
\title{Residual Plots for Linear and Generalized Linear Models}
\description{
Draws a plot or plots of residuals versus one or more term in a mean function and/or versus
fitted values.  For linear models curvature tests are computed for each of the plots
by adding a quadratic term to the regression function and testing the quadratic to be zero.  This is Tukey's test for nonadditivity when plotting against fitted values.
}
\usage{
### This is a generic function with only one required argument:

residualPlots (model, ...)

\S3method{residualPlots}{default}(model, terms= ~ . , 
            layout=NULL, ask, main="", 
            fitted=TRUE, AsIs=TRUE, plot=TRUE, tests=TRUE, groups, ...)

\S3method{residualPlots}{lm}(model, ...)

\S3method{residualPlots}{glm}(model, ...)

### residualPlots calls residualPlot, so these arguments can be
### used with either function

residualPlot(model, ...)

\S3method{residualPlot}{default}(model, variable = "fitted", type = "pearson",
     groups, plot = TRUE, linear = TRUE,     
     quadratic = if(missing(groups)) TRUE else FALSE, 
     smooth=FALSE, id=FALSE,
     col = carPalette()[1], col.quad = carPalette()[2], pch=1,
     xlab, ylab, lwd=1, grid=TRUE, key=!missing(groups), ...)

\S3method{residualPlot}{lm}(model, ...)

\S3method{residualPlot}{glm}(model, variable = "fitted", type = "pearson",
    plot = TRUE, quadratic = FALSE, smooth=TRUE, ...)
}
\arguments{
  \item{model}{
A regression object.
}
  \item{terms}{
A one-sided formula that specifies a subset of the terms that appear in the formula that defined the model. The default
\code{~ .} is to plot against all first-order terms.  Interactions are skipped.  For example, the
specification \code{terms = ~ . - X3} would plot against all terms
except for \code{X3}.  To get a plot against fitted values only, use the
arguments \code{terms = ~ 1}.  If a term like \code{log(X1)} is in the formula, then the corresponding
plot is obtained using \code{terms = ~ log(X1)}.
For polynomial terms generated using the \code{poly} function, the
plot is against the first-order variable, which may be centered and scaled
depending on how the arguments to the \code{poly} function. Plots against factors
are boxplots.  Plots against splines use the
result of \code{predict(model), type="terms")[, variable])} as the
horizontal axis.

A grouping variable can also be specified, so, for example
\code{terms= ~ .|type} would use the factor \code{type} to set a different
color and symbol for each level of \code{type}.  Any fits in the plots will
also be done separately for each level of group.  This can be useful of finding interactions
between the grouping variable and the term on the horizontal axis.  The grouping variable can
also be set using the argument \code{goups="type"}.
}
  \item{layout}{
If set to a value like \code{c(1, 1)} or \code{c(4, 3)}, the layout
of the graph will have this many rows and columns per "page" with a prompt given the in console to change to the next page.  If not set, the program
will select an appropriate layout.  If the number of graphs exceed nine, you
must select the layout yourself, or you will get a maximum of nine per page.
If \code{layout=NA}, the function does not set the layout and the user can
use the \code{par} function to control the layout, for example to have
plots from two models in the same graphics window.
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
Some terms in a model formula use the as-is or identity function, for example a term like \code{I(X^2)} to include a quadratic.  These terms will be skipped in plotting if \code{AsIs} is \code{FALSE}.  The default is \code{TRUE}.
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
Quoted variable name for the factor or regressor to be put on the horizontal axis, or
the default \code{"fitted"} to plot versus fitted values.
}
  \item{type}{
Specifies the type of residual to be plotted.  Any of \code{c("working", "response", "deviance", "pearson", "partial", "rstudent", "rstandard")} may be specified.  The default \code{type = "pearson"} is usually appropriate, since it is equal to ordinary residuals observed minus fit with ols, and correctly weighted residuals with wls or for a glm.  The last two options use the \code{\link{rstudent}} and \code{\link{rstandard}} functions and use studentized or standardized residuals.
}
  \item{groups}{A grouping indicator, provided as the quoted name of a grouping variable in the model, such as "type", or "ageGroup".  Points in different
groups will be plotted with different colors and symbols. If missing, no grouping is used.  In \code{residualPlots}, the grouping variable can also be set in the \code{terms}
argument, as described above.  The default is no grouping. 
}
  \item{linear}{If \code{TRUE}, adds a horizontal line at zero if no groups.
With groups, within level of groups display the ols regression of the residuals
as response and the horizontal axis as the regressor.
}
  \item{quadratic}{
if \code{TRUE}, fits the quadratic regression of the
vertical axis on the horizontal axis and displays a lack of fit test.  Default
is \code{TRUE} for \code{lm} and \code{FALSE} for \code{glm} or if \code{groups}
not missing.
}
  \item{smooth}{specifies the smoother to be used along with its arguments; if \code{FALSE}, which is the default except for
    generalized linear models, no smoother is shown;
    can be a list giving the smoother function and its named arguments; \code{TRUE} is equivalent to
    \code{list(smoother=loessLine, span=2/3, col=carPalette()[3])}, which is the default for a GLM.
    See \code{\link{ScatterplotSmoothers}} for the smoothers supplied by the
    \pkg{car} package and their arguments.}
  \item{id}{controls point identification; if \code{FALSE} (the default), no points are identified;
    can be a list of named arguments to the \code{\link{showLabels}} function;
    \code{TRUE} is equivalent to \code{list(method="r", n=2, cex=1, col=carPalette()[1], location="lr")},
    which identifies the 2 points with the largest absolute residuals.}
  \item{col}{
default color for points.  If groups is set, col can be a list at least
as long as the number of levels for groups giving the colors for each groups.
}
  \item{col.quad}{
default color for quadratic fit if groups is missing.  Ignored if groups are used.
}
  \item{pch}{plotting character.  The default is pch=1.  If groups are used,
pch can be set to a vector at least as long as the number of groups.
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
line width for the quadratic fit, if any.
}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
  \item{key}{Should a key be added to the plot?  Default is \code{!is.null(groups)}.
  }
}
\details{
\code{residualPlots} draws one or more residuals plots depending on the
value of the \code{terms} and \code{fitted} arguments.  If \code{terms = ~ .},
the default, then a plot is produced of residuals versus each first-order
term in the formula used to create the model.  A plot of
residuals versus fitted values is also included unless \code{fitted=FALSE}.  Setting \code{terms = ~1}
will provide only the plot against fitted values.

A table of curvature tests is displayed for linear models. For plots
against a term in the model formula, say \code{X1}, the test displayed is
the t-test for for \code{I(X1^2)} in the fit of \code{update, model, ~. + I(X1^2))}.
Econometricians call this a specification test.  For factors, the displayed
plot is a boxplot, no curvature test is computed, and grouping is ignored.
For fitted values in a linear model, the test is Tukey's one-degree-of-freedom test for
nonadditivity.  You can suppress the tests with the argument \code{tests=FALSE}.
If grouping is used curvature tests are not displayed.

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

Fox, J. and Weisberg, S. (2019)
\emph{An R Companion to Applied Regression}, Third Edition.  Sage.

Weisberg, S. (2014) \emph{Applied
Linear Regression}, Fourth Edition, Wiley, Chapter 8}
\author{Sanford Weisberg, \email{sandy@umn.edu}}

\seealso{See Also \code{\link{lm}}, \code{\link{identify}},
\code{\link{showLabels}} }


\examples{
m1 <- lm(prestige ~ income + type, data=Prestige)
residualPlots(m1)
residualPlots(m1, terms= ~ 1 | type) # plot vs. yhat grouping by type
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
