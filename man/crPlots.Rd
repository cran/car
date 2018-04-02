\name{crPlots}
\alias{crPlots}
\alias{crp}
\alias{crPlot}
\alias{crPlot.lm}

\title{Component+Residual (Partial Residual) Plots}
\description{
  These functions construct component+residual plots, also called
  partial-residual plots, for linear and generalized linear models.
}
\usage{
crPlots(model, terms = ~., layout = NULL, ask, main, 
    ...)

crp(...)

crPlot(model, ...)

\method{crPlot}{lm}(model, variable, id=FALSE,
    order=1, line=TRUE, smooth=TRUE, 
    col=carPalette()[1], col.lines=carPalette()[-1],
    xlab, ylab, pch=1, lwd=2, grid=TRUE, ...)
}

\arguments{
  \item{model}{model object produced by \code{lm} or \code{glm}.}
  \item{terms}{ A one-sided formula that specifies a subset of the regressors. 
    One component-plus-residual plot is drawn for each regressor.  The default
    \code{~.} is to plot against all numeric regressors.  For example, the
    specification \code{terms = ~ . - X3} would plot against all regressors
    except for \code{X3}, while \code{terms = ~ log(X4)} would give the plot for the 
    predictor X4 that is represented in the model by log(X4).  If this argument is a 
    quoted name of one of the predictors, the 
    component-plus-residual plot is drawn for that predictor only.}
  \item{layout}{ 
    If set to a value like \code{c(1, 1)} or \code{c(4, 3)}, the layout
    of the graph will have this many rows and columns.  If not set, the program
    will select an appropriate layout.  If the number of graphs exceed nine, you 
    must select the layout yourself, or you will get a maximum of nine per page.
    If \code{layout=NA}, the function does not set the layout and the user can
    use the \code{par} function to control the layout, for example to have 
    plots from two models in the same graphics window.}
  \item{ask}{If \code{TRUE}, ask the user before drawing the next plot; if \code{FALSE}, the default,  don't ask.
    This is relevant only if not all the graphs can be drawn in one window.}
  \item{main}{The title of the plot; if missing, one will be supplied.}
  \item{\dots}{\code{crPlots} passes these arguments to \code{crPlot}.  
    \code{crPlot} passes them to \code{plot}. }
  \item{variable}{A quoted string giving the name of a variable for the horizontal
    axis.}
  \item{id}{controls point identification; if \code{FALSE} (the default), no points are identified;
    can be a list of named arguments to the \code{\link{showLabels}} function;
    \code{TRUE} is equivalent to \code{list(method=list(abs(residuals(model, type="pearson")), "x"), n=2, 
    cex=1, col=carPalette()[1], location="lr")},
    which identifies the 2 points with the largest residuals and the 2 points with the most extreme
    horizontal (X) values.}
  \item{order}{order of polynomial regression performed for predictor to be plotted; default \code{1}.}
  \item{line}{\code{TRUE} to plot least-squares line.}
  \item{smooth}{specifies the smoother to be used along with its arguments; if \code{FALSE}, no smoother is shown;
    can be a list giving the smoother function and its named arguments; \code{TRUE}, the default, is equivalent to
    \code{list(smoother=loessLine)}. See \code{\link{ScatterplotSmoothers}} for the smoothers supplied by the
    \pkg{car} package and their arguments.}
  \item{col}{color for points; the default is the first entry
    in the current \pkg{car} palette (see \code{\link{carPalette}}
    and \code{\link[graphics]{par}}).}
  \item{col.lines}{a list of at least two colors.  The first color is used for the
    ls line and the second color is used for the fitted lowess line.  To use
    the same color for both, use, for example, \code{col.lines=c("red", "red")}}
  \item{xlab,ylab}{labels for the x and y axes, respectively.  If not set 
    appropriate labels are created by the function.}
  \item{pch}{plotting character for points; default is \code{1} 
    (a circle, see \code{\link[graphics]{par}}).}
  \item{lwd}{line width; default is \code{2} (see \code{\link[graphics]{par}}).}
  \item{grid}{If TRUE, the default, a light-gray background grid is put on the
    graph.}
}

\details{
  The function intended for direct use is \code{crPlots}, for which \code{crp}
  is an abbreviation. 
  
  The model cannot contain interactions, but can contain factors. 
  Parallel boxplots of the partial residuals are drawn for the levels
  of a factor.
}
\value{
  \code{NULL}. These functions are used for their side effect of producing
  plots.
}
\references{ 
  Cook, R. D. and Weisberg, S. (1999)
  \emph{Applied Regression, Including Computing and Graphics.} Wiley.

  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{ceresPlots}}, \code{\link{avPlots}}}

\examples{
crPlots(m<-lm(prestige ~ income + education, data=Prestige)) 

crPlots(m, terms=~ . - education) # get only one plot

crPlots(lm(prestige ~ log2(income) + education + poly(women,2), data=Prestige))

crPlots(glm(partic != "not.work" ~ hincome + children, 
  data=Womenlf, family=binomial), smooth=list(span=0.75))
}

\keyword{hplot}
\keyword{regression}
