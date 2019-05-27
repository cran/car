\name{mcPlots}
\alias{mcPlots}
\alias{mcPlot}
\alias{mcPlot.lm}


\title{Draw Linear Model Marginal and Conditional Plots in Parallel or Overlaid}
\description{
the \code{mcPlot} function draws two graphs or overlays the two graphs.  For a response Y and a regressor X, the first plot is the \emph{m}arginal plot of Y versus X with both variables centered, visualizing the conditional distribution of Y given X ignoring all other regressors.  The second plot is an added-variable for X after all other regressors, visualizing the \emph{c}onditional distribution of Y given X after adjusting for all other predictors.  The added variable plot by default is drawn using the same xlim and ylim as the centered marginal plot to emphasize that conditioning removes variation in both the regressor and the response.The plot is primarily intended as a pedagogical tool for understanding coefficients in first-order models.
}

\usage{
mcPlots(model, terms=~., layout=NULL, ask, overlaid=TRUE, ...)

mcPlot(model, ...)

\method{mcPlot}{lm}(model, variable, id=FALSE,
    col.marginal=carPalette()[2], col.conditional=carPalette()[3],
    col.arrows="gray", pch = c(16,1), lwd = 2, grid=TRUE,
    ellipse=FALSE, overlaid=TRUE, new=TRUE, title=TRUE, ...)
}
\arguments{
  \item{model}{model object produced by \code{lm}.
}
  \item{terms}{ A one-sided formula that specifies a subset of the predictors.
One added-variable plot is drawn for each regressor and for each basis vector used to define a factor. For example, the
specification \code{terms = ~ . - X3} would plot against all terms
except for \code{X3}. If this argument is a quoted name of one of the regressors or factors, the
added-variable plot is drawn for that regressor or factor only.  Unlike other car functions, the formula should include the names of regressors, not predictors.  That is, if \code{log(X4)} is used to represent a predictor \code{X4}, the formula should specify \code{terms = ~ log(X4)}.
}
\item{variable}{A quoted string giving the name of a numeric predictor in the model matrix for the horizontal
  axis.  To plot against a factor, you need to specify the full name of one of the indicator variables that define the factor.  For example, for a factor called \code{type} with levels \code{A}, \code{B} and {C}, using the usual drop-first level parameterization of the factor, the regressors for \code{type} would be \code{typeB} or \code{typeC}.  Similarly, to plot against the regressor \code{log(X4)}, you must specify \code{"log((X4)"}, not \code{"X4"}.
}
  \item{layout}{
If set to a value like \code{c(1, 2)} or \code{c(6, 2)}, the layout
of the graph will have this many rows and columns.  If not set, behavior depends on the value of the \code{overlaid} argument; see the details
}
\item{ask}{If \code{TRUE}, ask the user before drawing the next plot; if \code{FALSE} don't ask.
}
\item{\dots}{\code{mcPlots} passes these arguments to \code{mcPlot}.
  \code{mcPlot} passes arguments to \code{plot}.
}
\item{id}{controls point identification; if \code{FALSE} (the default), no points are identified;
    can be a list of named arguments to the \code{\link{showLabels}} function;
    \code{TRUE} is equivalent to \code{list(method=list(abs(residuals(model, type="pearson")), "x"), n=2, cex=1, col=carPalette()[1], location="lr")},
    which identifies the 2 points with the largest residuals and the 2 points with the most extreme
    horizontal (X) values.}
  \item{overlaid}{If TRUE, the default, overlay the marginal and conditional plots on the same graph; otherwise
  plot them side-by-side.  See the details below}
  \item{col.marginal, col.conditional}{colors for points, lines,
  ellipses in the marginal and conditional plots, respectively.  The defaults are determined by the \code{\link{carPalette}} function.}
  \item{col.arrows}{color for the arrows with \code{overlaid=TRUE}}
  \item{pch}{Plotting character for marginal and conditional plots, respectively.}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{grid}{If \code{TRUE}, the default, a light-gray background grid is put on the graph.}
  \item{ellipse}{Arguments to pass to the \code{\link{dataEllipse}} function, in the form of a list
    with named elements; e.g., \code{ellipse.args=list(robust=TRUE))} will cause the ellipse to be plotted using
    a robust covariance-matrix. if \code{FALSE}, the default, no ellipse is plotted. \code{TRUE} is equivalent to
    \code{ellipse=list(levels=0.5)}, which plots a bivariate-normal 50 percent concentration ellipse.}
  \item{new}{if \code{TRUE}, the default, the plot window is reset when \code{overlaid=FALSE} using \code{par{mfrow=c(1, 2)}}.  If \code{FALSE}, the layout of the plot window is not reset.
  Users will ordinarily ignore this argument.}
  \item{title}{If TRUE, the default, the standard main argument in plot is used to add a standard title to each plot.  If FALSE no title is used.}
}

\details{
With an \code{lm} object, suppose the response is Y, X is a numeric regressor of interest, and Z is all the remaining predictors, possibly including interactions and factors.  This function produces two graphs.  The first graph is the marginal plot of Y versus X, with each variable centered around its mean.  The second conditional plot is the added-variable plot of e(Y|Z) versus e(X|Z) where e(a|b) means the Pearson residuals from the regression of a on b.  If \code{overlaid=TRUE}, these two plots are overlaid in one graph, with the points in different colors.  In addition, each point in the marginal plot is joined to its value in the conditional plot by an arrow.  Least squares regression lines fit to the marginal and conditional graphs are also shown; data ellipsoids can also be added.  If \code{overlaid=FALSE}, then the two graphs are shown in side-by-side plots as long as the second argument to \code{layout} is equal to \code{2}, or \code{layout} is set by the function. The arrows are omitted if the graphs are not overlaid.

These graphs are primarily for teaching, as the marginal plot shows the relationship between Y and X ignoring Z, while the conditional is the relationship between Y and X given X.  By keeping the scales the same in both graphs the effect of conditioning on both X and Y can be visualized.

This function is intended for first-order models with numeric predictors only.  For a factor,  one (pair) of mcPlots will be produced for each of the dummy variables in the basis for the factor, and the resulting plots are not generally meaningful because they depend on parameterization.  If the mean function includes interactions, then mcPlots for main effects may violate the hierarchy principle, and may also be of little interest.  mcPlots for interactions of numerical predictors, however, can be useful.

These graphs are closely related to the ARES plots proposed by Cook and Weisberg (1989).
This plot would benefit from animation.
}

\value{
  These functions are used for their side effect of producing
  plots.
}

\references{
  Cook, R. D. and Weisberg, S. (1989)
  \emph{Regression diagnostics with dynamic graphics,} Technometrics, 31, 277.

  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.

  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.

  Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley.
}

\author{John Fox \email{jfox@mcmaster.ca}, Sanford Weisberg \email{sandy@umn.edu}}

\seealso{\code{\link{avPlots}}, \code{\link{residualPlots}}, \code{\link{crPlots}}, \code{\link{ceresPlots}},
\code{\link{dataEllipse}}
}

\examples{
m1 <- lm(partic ~ tfr + menwage + womwage + debt + parttime, data = Bfox)
mcPlot(m1, "womwage")
mcPlot(m1, "womwage", overlaid=FALSE, ellipse=TRUE)
}

\keyword{hplot}
\keyword{regression}
