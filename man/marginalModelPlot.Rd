\name{mmps}
\alias{mmps}
\alias{mmp}
\alias{mmp.lm}
\alias{mmp.glm}
\alias{mmp.default}
\alias{marginalModelPlot}
\alias{marginalModelPlots}

\title{Marginal Model Plotting}
\description{
For a regression object, draw a plot of the response on the vertical axis versus
a linear combination \eqn{u} of regressors in the mean function on the horizontal
axis.  Added to the plot are a smooth for the graph, along with
a smooth from the plot of the fitted values on \eqn{u}.  \code{mmps} is an alias
for \code{marginalModelPlots}, and \code{mmp} is an alias for \code{marginalModelPlot}.
}
\usage{
marginalModelPlots(...)

mmps(model, terms= ~ ., fitted=TRUE, layout=NULL, ask,
        main, groups, key=TRUE, ...)

marginalModelPlot(...)

mmp(model, ...)

\method{mmp}{lm}(model, variable, sd = FALSE,
    xlab = deparse(substitute(variable)), 
    smooth=TRUE, key=TRUE, pch, groups=NULL, ...)

\method{mmp}{default}(model, variable, sd = FALSE,
    xlab = deparse(substitute(variable)), ylab, smooth=TRUE,
    key=TRUE, pch, groups=NULL,
    col.line = carPalette()[c(2, 8)], col=carPalette()[1],
    id=FALSE, grid=TRUE, ...)

\method{mmp}{glm}(model, variable, sd = FALSE,
    xlab = deparse(substitute(variable)), ylab,
    smooth=TRUE, key=TRUE, pch, groups=NULL,
    col.line = carPalette()[c(2, 8)], col=carPalette()[1],
    id=FALSE, grid=TRUE, ...)
}

\arguments{
\item{model}{A regression object, usually of class either \code{lm} or \code{glm},
   for which there is a \code{predict} method defined. }
\item{terms}{A one-sided formula.  A marginal model plot will be drawn for
  each term on the right-side of this formula that is not a factor.  The
  default is \code{~ .}, which specifies that all the terms in
  \code{formula(object)} will be used.  If a conditioning argument is given,
  eg \code{terms = ~. | a}, then separate colors and smoothers are used for
  each unique non-missing value of \code{a}.  See examples below.}
\item{fitted}{If \code{TRUE}, the default, then a marginal model plot in the direction
  of the fitted values for a linear model or the linear predictor of a 
  generalized linear model will be drawn.}
\item{layout}{
  If set to a value like \code{c(1, 1)} or \code{c(4, 3)}, the layout
of the graph will have this many rows and columns.  If not set, the program
will select an appropriate layout.  If the number of graphs exceed nine, you
must select the layout yourself, or you will get a maximum of nine per page.
If \code{layout=NA}, the function does not set the layout and the user can
use the \code{par} function to control the layout, for example to have
plots from two models in the same graphics window.
}
\item{ask}{If \code{TRUE}, ask before clearing the graph window to draw more plots.}
        \item{main}{
Main title for the array of plots.  Use \code{main=""} to suppress the title;
if missing, a title will be supplied.
}
\item{\dots}{
Additional arguments passed from \code{mmps} to \code{mmp} and
then to \code{plot}.  Users should generally use \code{mmps}, or equivalently
\code{marginalModelPlots}.
}
\item{variable}{ The quantity to be plotted on the horizontal axis.  If this argument
  is missing, the horizontal variable is the linear predictor, returned by 
  \code{predict(object)} for models of class
  \code{lm}, with default label \code{"Fitted values"},  or returned by 
  \code{predict(object, type="link")} for models of class \code{glm}, with default 
  label \code{"Linear predictor"}. It can be any other
  vector of length equal to the number of observations in the object. Thus the
  \code{mmp} function can be used to get a marginal model plot versus any
  regressor or predictor while the \code{mmps} function can be used only to get
  marginal model plots for the first-order regressors in the formula.  In
  particular, terms defined by a spline basis are skipped by \code{mmps}, but
  you can use \code{mmp} to get the plot for the variable used to define
  the splines.}
\item{sd}{ If \code{TRUE}, display sd smooths.  For a binomial regression with all
sample sizes equal to one, this argument is ignored as the SD bounds don't
make any sense. }
\item{xlab}{ label for horizontal axis.}
\item{ylab}{ label for vertical axis, defaults to name of response.}
\item{smooth}{specifies the smoother to be used along with its arguments; if \code{FALSE}, no smoother is shown;
can be a list giving the smoother function and its named arguments; \code{TRUE}, the default, is equivalent to
\code{list(smoother=loessLine, span=2/3)} for linear models and \code{list(smoother=gamLine, k=3)} for generalized linear models.
See \code{\link{ScatterplotSmoothers}} for the smoothers supplied by the
\pkg{car} package and their arguments; the \code{spread} argument is not supported for marginal model plots.}
\item{groups}{The name of a vector that specifies a grouping variable for
separate colors/smoothers.  This can also be specified as a conditioning
argument on the \code{terms} argument.}
\item{key}{If \code{TRUE}, include a key at the top of the plot, if \code{FALSE} omit the
key.  If grouping is present, the key is only printed for the upper-left plot.}
\item{id}{controls point identification; if \code{FALSE} (the default), no points are identified;
can be a list of named arguments to the \code{\link{showLabels}} function;
\code{TRUE} is equivalent to \code{list(method="y", n=2, cex=1, col=carPalette()[1], location="lr")},
which identifies the 2 points with the most unusual response (Y) values.}
\item{pch}{plotting character to use if no grouping is present.}
\item{col.line}{ colors for data and model smooth, respectively.  The default is to use \code{\link{carPalette}}, \code{carPalette()[c(2, 8)]},  blue and red. }
\item{col}{color(s) for the plotted points.}
\item{grid}{If TRUE, the default, a light-gray background grid is put on the
graph}
}
\details{
\code{mmp} and \code{marginalModelPlot} draw one marginal model plot against
whatever is specified as the horizontal axis.
\code{mmps} and \code{marginalModelPlots} draws marginal model plots
versus each of the terms in the \code{terms} argument and versus fitted values.
\code{mmps} skips factors and  interactions if they are specified in the
\code{terms} argument.  Terms based on polynomials or on splines (or
potentially any term that is represented by a matrix of regressors) will
be used to form a marginal model plot by returning a linear combination of the
terms.  For example, if you specify \code{terms = ~ X1 + poly(X2, 3)} and
\code{poly(X2, 3)} was part of the original model formula, the horizontal
axis of the marginal model plot for \code{X2} will be the value of
\code{predict(model, type="terms")[, "poly(X2, 3)"])}.  If the \code{predict}
method for the model you are using doesn't support \code{type="terms"},
then the polynomial/spline term is skipped.  Adding a conditioning variable,
e.g., \code{terms = ~ a + b | c}, will produce marginal model plots for \code{a}
and \code{b} with different colors and smoothers for each unique non-missing
value of \code{c}.

For linear models, the default smoother is loess.
For generalized linear models, the default smoother uses \code{gamLine}, fitting
a generalized additive model with the same family, link and weights as the fit of the
model. SD smooths are not computed for for generalized linear models.

For generalized linear models the default number of elements in the spline basis is
\code{k=3}; this is done to allow fitting for predictors with just a few support
points.  If you have many support points you may wish to set \code{k} to a higher
number, or \code{k=-1} for the default used by \code{\link{gam}}.
}

\value{
Used for its side effect of producing plots.
}

\seealso{\code{\link{ScatterplotSmoothers}}, \code{\link{plot}} }

\references{

Cook, R. D., & Weisberg, S. (1997). Graphics for assessing the adequacy of regression models.
\emph{Journal of the American Statistical Association}, 92(438), 490-499.

Fox, J. and Weisberg, S. (2019)
\emph{An R Companion to Applied Regression}, Third Edition.  Sage.

Weisberg, S. (2005) \emph{Applied
Linear Regression}, Third Edition, Wiley, Section 8.4.
}

\author{Sanford Weisberg, \email{sandy@umn.edu}}

\examples{
c1 <- lm(infantMortality ~ ppgdp, UN)
mmps(c1)
c2 <- update(c1, ~ log(ppgdp))
mmps(c2)
# include SD lines
p1 <- lm(prestige ~ income + education, Prestige)
mmps(p1, sd=TRUE)
# condition on type:
mmps(p1, ~. | type)
# logisitic regression example
# smoothers return warning messages.
# fit a separate smoother and color for each type of occupation.
m1 <- glm(lfp ~ ., family=binomial, data=Mroz)
mmps(m1)
}
\keyword{hplot}
\keyword{regression}
