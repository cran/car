\name{ScatterplotSmoothers}
\alias{ScatterplotSmoothers}
\alias{gamLine}
\alias{quantregLine}
\alias{loessLine}

\title{
Smoothers to Draw Lines on Scatterplots
}

\description{
These smoothers are used to draw nonparametric-regression lines on scatterplots produced by
the \code{\link{scatterplot}}, \code{\link{scatterplotMatrix}} and other \code{car} functions. 
The functions aren't meant to
be called directly by the user, although the user can supply options via the \code{smoother.args} argument,
the contents of which vary by the smoother (see \emph{Details} below). The \code{gamLine} smoother uses the
\code{\link[mgcv]{gam}} function in the \pkg{mgcv} package, the \code{loessLine} smoother uses the
\code{\link[stats]{loess}} function in the \pkg{stats} package, and the \code{quantregLine} smoother uses the
\code{\link[quantreg]{rqss}} function in the \pkg{quantreg} package.
}

\usage{
gamLine(x, y, col, log.x, log.y, spread=FALSE, smoother.args, draw=TRUE, offset=0)

loessLine(x, y, col, log.x, log.y, spread=FALSE, smoother.args, draw=TRUE, offset=0)

quantregLine(x, y, col, log.x, log.y, spread=FALSE, smoother.args, draw=TRUE, offset=0)
}

\arguments{
  \item{x}{$x$ coordinates of points.}
  \item{y}{$y$ coordinates of points.}
  \item{col}{line color.}
  \item{log.x}{\code{TRUE} if the $x$-axis is logged.}
  \item{log.y}{\code{TRUE} if the $y$-axis is logged.}
  \item{spread}{the default is to plot only an estimated mean or median function.  If this
  argument is TRUE, then a measure of spread is also plotted.}
  \item{smoother.args}{additional options accepted by the smoother, in the form of a list of
    named values (see \emph{Details} below).}
  \item{draw}{if TRUE, the default, draw the smoother on the currently active graph.
    If FALSE, return a list with coordinates \code{x} and \code{y} for the points 
    that make up the smooth and if requested \code{x.pos, y.pos, x.neg, y.neg} for the
    spread smooths.}
  \item{offset}{For use when \code{spread=TRUE}, the vertical axis is \code{sqrt(offset^2 + variance smooth)}.}
}

\details{
The function \code{loessLine} is a re-implementation of the \code{loess} smoother
that was used in \code{car} prior to September 2012.  The only enhancement is the ability to
set more arguments through the \code{smoother.args} argument.

The function \code{gamLine} is new and more general than the \code{loess} fitting
because it allows fitting a generalized additive model using splines.  You can specify an error
distribution and link function.  

The function \code{quantregLine} fits an additive model using splines with estimation
based on L1 regression and quantile regression if you ask for the spread.  It is
likely to be more robust than the other smoothers.

The argument \code{smoother.args} is a list of named elements used to pass 
additional arguments to the smoother.

For \code{loessLine} the default value is
\code{smoother.args=list(lty=1, lwd=2, lty.spread=2, lwd.spread=1, span=0.5,
degree=2, family="symmetric", iterations=4)}.  
The arguments \code{lty} and \code{lwd} are the type and width 
respectively of the mean or median smooth, \code{smooth.lty}
and \code{smooth.lwd} are the type and color of the spread smooths if requested.
The arguments \code{span}, \code{degree} and \code{family} are
passed to the \code{\link{loess}} function, \code{iterations=0} by default
specifies no robustness iterations.

For \code{gamLine} the default is
\code{smoother.args=list(lty=1, lwd=2, lty.spread=2, lwd.spread=1,
k=-1, bs="tp", family="gaussian", link=NULL, weights=NULL)}
The first four arguments are as for \code{loessLine}.  The next two 
arguments are passed to the \code{gam} function to control the smoothing:
\code{k=-1} allows \code{gam} to choose the number of splines in the basis
function; \code{bs="tp"} provides the type of spline basis to be used with \code{"tp"} 
for the default thin-plate splines.  The last three arguments allow providing
a family, link and weights as in generalized linear models.  See examples
below.  The spread argument is ignored unless family="gaussian" and link=NULL.

For \code{quantregLine} the default is
\code{smoother.args=list(lty=1, lwd=2, lty.spread=2, lwd.spread=1,
lambda=IQR(x))}.  The first four
arguments are as for \code{loessLine}.  The last argument is passed to the
\code{\link[quantreg]{qss}} function in \code{quantreg}.  It is a smoothing
parameter, here a robust estimate of the scale of the horizontal axis variable.
This is an arbitrary choice, and may not work well in all circumstances.
}
    
\author{John Fox \email{jfox@mcmaster.ca} and Sanford Weisberg \email{sandy@umn.edu}.}

\seealso{\code{\link{scatterplot}}, \code{\link{scatterplotMatrix}}, \code{\link[mgcv]{gam}},
    \code{\link[stats]{loess}}, and \code{\link[quantreg]{rqss}}.}
    
\examples{
scatterplot(prestige ~ income, data=Prestige)
scatterplot(prestige ~ income, data=Prestige, smoother=gamLine)
scatterplot(prestige ~ income, data=Prestige, smoother=quantregLine)

scatterplot(prestige ~ income | type, data=Prestige)
scatterplot(prestige ~ income | type, data=Prestige, smoother=gamLine)
scatterplot(prestige ~ income | type, data=Prestige, smoother=quantregLine)
scatterplot(prestige ~ income | type, data=Prestige, smoother=NULL)

scatterplot(prestige ~ income | type, data=Prestige, spread=TRUE)
scatterplot(prestige ~ income | type, data=Prestige, smoother=gamLine, spread=TRUE)
scatterplot(prestige ~ income | type, data=Prestige, smoother=quantregLine, spread=TRUE)

scatterplot(weight ~ repwt | sex, spread=TRUE, data=Davis, smoother=loessLine)
scatterplot(weight ~ repwt | sex, spread=TRUE, data=Davis, smoother=gamLine) # messes up
scatterplot(weight ~ repwt | sex, spread=TRUE, data=Davis, smoother=quantregLine) #  robust
set.seed(12345)
w <- 1 + rpois(100, 5)
x <- rnorm(100)
p <- 1/(1 + exp(-(x + 0.5*x^2)))
y <- rbinom(100, w, p)
scatterplot(y/w ~ x, smoother=gamLine, 
       smoother.args=list(family="binomial", weights=w))
scatterplot(y/w ~ x, smoother=gamLine, 
       smoother.args=list(family=binomial, link="probit", weights=w))
scatterplot(y/w ~ x, smoother=gamLine, 
       smoother.args=list(family=binomial, link="probit", weights=w)) 
scatterplot(y/w ~ x, smoother=loessLine, reg=FALSE) 

y <- rbinom(100, 1, p)
scatterplot(y ~ x, smoother=gamLine, smoother.args=list(family=binomial)) 
}
\keyword{hplot}
