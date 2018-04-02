\name{dfbetaPlots}
\alias{dfbetaPlots}
\alias{dfbetasPlots}
\alias{dfbetaPlots.lm}
\alias{dfbetasPlots.lm}


\title{dfbeta and dfbetas Index Plots}
\description{
  These functions display index plots of dfbeta (effect on coefficients of deleting
  each observation in turn) and dfbetas (effect on coefficients of deleting each observation
  in turn, standardized by a deleted estimate of the coefficient standard error). In the
  plot of dfbeta, horizontal lines are drawn at 0 and +/- one standard error; in the plot
  of dfbetas, horizontal lines are drawn and 0 and +/- 1.
}

\usage{
dfbetaPlots(model, ...)

dfbetasPlots(model, ...)

\method{dfbetaPlots}{lm}(model, terms= ~ ., intercept=FALSE, layout=NULL, ask, 
	    main, xlab, ylab, labels=rownames(dfbeta), 
        id.method="y",  
        id.n=if(id.method[1]=="identify") Inf else 0, id.cex=1,
        id.col=carPalette()[1], id.location="lr", col=carPalette()[1], grid=TRUE, ...)

\method{dfbetasPlots}{lm}(model, terms=~., intercept=FALSE, layout=NULL, ask, 
	    main, xlab, ylab,
        labels=rownames(dfbeta), id.method="y",
        id.n=if(id.method[1]=="identify") Inf else 0, id.cex=1, 
        id.col=carPalette()[1], id.location="lr", col=carPalette()[1], grid=TRUE, ...)
}
	
\arguments{
  \item{model}{model object produced by \code{lm} or \code{glm}.
}
  \item{terms}{ A one-sided formula that specifies a subset of the terms in the model. 
One dfbeta or dfbetas plot is drawn for each regressor.  The default
\code{~.} is to plot against all terms in the model with the exception of an intercept.  For example, the
specification \code{terms = ~.-X3} would plot against all terms
except for \code{X3}.  If this argument is a quoted name of one of the terms, the 
index plot is drawn for that term only.
}
\item{intercept}{Include the intercept in the plots; default is \code{FALSE}.}
  \item{layout}{ 
If set to a value like \code{c(1, 1)} or \code{c(4, 3)}, the layout
of the graph will have this many rows and columns.  If not set, the program
will select an appropriate layout.  If the number of graphs exceed nine, you 
must select the layout yourself, or you will get a maximum of nine per page.
If \code{layout=NA}, the function does not set the layout and the user can
use the \code{par} function to control the layout, for example to have 
plots from two models in the same graphics window.
}
\item{main}{The title of the graph; if missing, one will be supplied.
}
\item{xlab}{Horizontal axis label; defaults to \code{"Index"}.}
\item{ylab}{Vertical axis label; defaults to coefficient name.}
\item{ask}{If \code{TRUE}, ask the user before drawing the next plot; if \code{FALSE}, the default, don't ask.
}
\item{\dots}{optional additional arguments to be passed to \code{\link{plot}}, \code{\link{points}}, and \code{\link{showLabels}}}. 
  \item{id.method, labels, id.n, id.cex, id.col, id.location}{Arguments for the labelling of
    points.  The default is \code{id.n=0} for labeling no points.  See
    \code{\link{showLabels}} for details of these arguments.
}
\item{col}{color for points; defaults to the first entry in the color \code{\link{carPalette}}.}
   \item{grid}{If \code{TRUE}, the default, a light-gray background grid is put on the
graph}
}


\value{
  \code{NULL}. These functions are used for their side effect: producing
  plots.
}

\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link{dfbeta}} ,\code{\link{dfbetas}}}

\examples{
dfbetaPlots(lm(prestige ~ income + education + type, data=Duncan))

dfbetasPlots(glm(partic != "not.work" ~ hincome + children, 
  data=Womenlf, family=binomial))
}

\keyword{hplot}
\keyword{regression}
