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
the \code{\link{scatterplot}}, \code{\link{scatterplotMatrix}} and several other \pkg{car} functions.
The functions are not meant to
be called directly by the user, although the user can supply options via the \code{smoother.args} argument,
the contents of which vary by the smoother (see \emph{Details} below). The \code{gamLine} smoother uses the
\code{\link[mgcv]{gam}} function in the \pkg{mgcv} package, the \code{loessLine} smoother uses the
\code{\link[stats]{loess}} function in the \pkg{stats} package, and the \code{quantregLine} smoother uses the
\code{\link[quantreg]{rqss}} function in the \pkg{quantreg} package.
}

\usage{
gamLine(x, y, col=carPalette()[1], log.x=FALSE, log.y=FALSE, var=FALSE, spread=var,
    smoother.args=NULL, draw=TRUE, offset=0)

loessLine(x, y, col=carPalette()[1], log.x=FALSE, log.y=FALSE, var=FALSE, spread=var,
    smoother.args=NULL, draw=TRUE, offset=0)

quantregLine(x, y, col=carPalette()[1], log.x=FALSE, log.y=FALSE, var=FALSE, spread=var,
    smoother.args=NULL, draw=TRUE, offset=0)
}

\arguments{
  \item{x}{horizontal coordinates of points.}
  \item{y}{vertical coordinates of points.}
  \item{col}{line color.}
  \item{log.x}{should be set to \code{TRUE} (default is \code{FALSE}) if the horizontal axis is logged.}
  \item{log.y}{should be set to \code{TRUE} (default is \code{FALSE}) if the vertical axis is logged.}
  \item{spread, var}{the default is to plot only an estimated mean or median function.  If either of
  these arguments is TRUE, then a measure of variability is also plotted.}
  \item{smoother.args}{additional options accepted by the smoother, in the form of a list of
    named values (see \emph{Details} below).}
  \item{draw}{if TRUE, the default, draw the smoother on the currently active graph.
    If FALSE, return a list with coordinates \code{x} and \code{y} for the points
    that make up the smooth and if requested \code{x.pos, y.pos, x.neg, y.neg} for the
    spread smooths.}
  \item{offset}{For use when \code{spread=TRUE}, the vertical axis is \code{sqrt(offset^2 + variance smooth)}.}
}

\details{
The function \code{loessLine} is a re-implementation of the \pkg{loess} smoother
that was used in \pkg{car} prior to September 2012.  The main enhancement is the ability to
set more arguments through the \code{smoother.args} argument.

The function \code{gamLine} is more general than the \code{loess} fitting
because it allows fitting a generalized additive model using splines.  You can specify an error
distribution and link function.

The function \code{quantregLine} fits an additive model using splines with estimation
based on L1 regression for the median and quantile regression if you ask for the spread.  It is
likely to be more robust than the other smoothers.

The argument \code{smoother.args} is a list of named elements used to pass
additional arguments to the smoother.  As of November, 2016, the smoother is evaluated at an equally spaced grid of 50 points in the range of the horizontal variable.  With any of the smoothers you can change to say 100 evaluation points by using the argument \code{smoother.args=list(evaluation=100)}.  As of version 3.0-1, the arguments in \code{smoother.args} \code{col.var}, \code{lty.var} and \code{lwd.var} are equivalent to \code{col.spread}, \code{lty.spread} and \code{lwd.spread}, respectively.

For \code{loessLine} the default value is
\code{smoother.args=list(lty.smooth=1, lwd.smooth=2, lty.spread=4, lwd.spread=2, span=2/3 (prior to 11/2016, span was 1/2),
degree=1, family="symmetric", iterations=4)}.
The arguments \code{lty.smooth} and \code{lwd.smooth} are the type and width
respectively of the mean or median smooth, \code{lty.spread}
and \code{lwd.spred} are the type and color of the spread smooths if requested.
The arguments \code{span}, \code{degree} and \code{family} are
passed to the \code{\link{loess}} function, \code{iterations=4} robustness iterations.

For \code{gamLine} the default is
\code{smoother.args=list(lty.smooth=1, lwd.smooth=2, lty.spread=4, lwd.spread=2,
k=-1, bs="tp", family="gaussian", link=NULL, weights=NULL)}
The first four arguments are as for \code{loessLine}.  The next two
arguments are passed to the \code{gam} function to control the smoothing:
\code{k=-1} allows \code{gam} to choose the number of splines in the basis
function; \code{bs="tp"} provides the type of spline basis to be used with \code{"tp"}
for the default thin-plate splines.  The last three arguments allow providing
a family, link and weights as in generalized linear models.  See examples
below.  The spread argument is ignored unless family="gaussian" and link=NULL.

For \code{quantregLine} the default is
\code{smoother.args=list(lty.smooth=1, lwd.smooth=2, lty.spread=4, lwd.spread=2,
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
scatterplot(prestige ~ income, data=Prestige, smooth=list(smoother=gamLine))
scatterplot(prestige ~ income, data=Prestige,
    smooth=list(smoother=quantregLine))

scatterplot(prestige ~ income | type, data=Prestige)
scatterplot(prestige ~ income | type, data=Prestige,
    smooth=list(smoother=gamLine))
scatterplot(prestige ~ income | type, data=Prestige,
    smooth=list(smoother=quantregLine))
scatterplot(prestige ~ income | type, data=Prestige, smooth=FALSE)

scatterplot(prestige ~ income | type, data=Prestige, spread=TRUE)
scatterplot(prestige ~ income | type, data=Prestige,
    smooth=list(smoother=gamLine), spread=TRUE)
scatterplot(prestige ~ income | type, data=Prestige,
    smooth=list(smoother=quantregLine), spread=TRUE)

scatterplot(weight ~ repwt | sex, spread=TRUE, data=Davis,
    smooth=list(smoother=loessLine))
scatterplot(weight ~ repwt | sex, spread=TRUE, data=Davis,
    smooth=list(smoother=gamLine)) # messes up
scatterplot(weight ~ repwt | sex, spread=TRUE, data=Davis,
    smooth=list(smoother=quantregLine)) #  robust

set.seed(12345)
w <- 1 + rpois(100, 5)
x <- rnorm(100)
p <- 1/(1 + exp(-(x + 0.5*x^2)))
y <- rbinom(100, w, p)
scatterplot(y/w ~ x, smooth=list(smoother=gamLine, family="binomial",
    weights=w))
scatterplot(y/w ~ x, smooth=list(smoother=gamLine, family=binomial,
    link="probit", weights=w))
scatterplot(y/w ~ x, smooth=list(smoother=loessLine), reg=FALSE)

y <- rbinom(100, 1, p)
scatterplot(y ~ x, smooth=list(smoother=gamLine, family=binomial))
}
\keyword{hplot}
