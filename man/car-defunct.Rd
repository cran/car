\name{car-defunct}
\alias{car-defunct}
\alias{av.plot}
\alias{av.plots}
\alias{box.cox}
\alias{bc}
\alias{box.cox.powers}
\alias{box.cox.var}
\alias{box.tidwell}
\alias{cookd}
\alias{confidence.ellipse}
\alias{ceres.plot}
\alias{ceres.plots}
\alias{cr.plot}
\alias{cr.plots}
\alias{data.ellipse}
\alias{durbin.watson}
\alias{levene.test}
\alias{leverage.plot}
\alias{leverage.plots}
\alias{linear.hypothesis}
\alias{outlier.test}
\alias{ncv.test}
\alias{qq.plot}
%\alias{scatterplot.matrix}
\alias{spread.level.plot}
\alias{skewPower}

\title{Defunct Functions in the car Package}

\description{
  These functions are were deprecated in 2009 and are now defunct.
}

\usage{
av.plot(...)
av.plots(...)
box.cox(...)
bc(...)
box.cox.powers(...)
box.cox.var(...)
box.tidwell(...)
cookd(...)
confidence.ellipse(...)
ceres.plot(...)
ceres.plots(...)
cr.plot(...)
cr.plots(...)
data.ellipse(...)
durbin.watson(...)
levene.test(...)
leverage.plot(...)
leverage.plots(...)
linear.hypothesis(...)
ncv.test(...)
outlier.test(...)
qq.plot(...)
skewPower(...)
spread.level.plot(...)
}

\arguments{
  \item{\dots}{pass arguments down.}
}

\details{
\code{av.plot} and \code{av.plots} are replaced by \code{\link{avPlot}} and \code{\link{avPlots}} functions.

\code{box.cox} and \code{bc} are now replaced by \code{\link{bcPower}}.

\code{box.cox.powers} is replaced by \code{\link{powerTransform}}.

\code{box.cox.var} is replaced by \code{\link{boxCoxVariable}}.

\code{box.tidwell} is replaced by \code{\link{boxTidwell}}.

\code{cookd} is replaced by \code{\link[stats:influence.measures]{cooks.distance}} in the \pkg{stats} package.

\code{confidence.ellipse} is replaced by \code{\link{confidenceEllipse}}.

\code{ceres.plot} and \code{ceres.plots} are now replaced by the \code{\link{ceresPlot}} and \code{\link{ceresPlots}} 
functions.

\code{cr.plot} and \code{cr.plots} are now replaced by the \code{\link{crPlot}} and \code{\link{crPlots}} functions.

\code{data.ellipse} is replaced by \code{\link{dataEllipse}}.

\code{durbin.watson} is replaced by \code{\link{durbinWatsonTest}}.

\code{levene.test} is replaced by \code{\link{leveneTest}} function.

\code{leverage.plot} and \code{leverage.plots} are now replaced by the \code{\link{leveragePlot}} and \code{\link{leveragePlots}} 
functions.

\code{linear.hypothesis} is replaced by the \code{\link{linearHypothesis}} function.

\code{ncv.test} is replaced by \code{\link{ncvTest}}.

\code{outlier.test} is replaced by \code{\link{outlierTest}}.

\code{qq.plot} is replaced by \code{\link{qqPlot}}.

\code{skewPower} is replaced by \code{\link{bcnPower}}.

\code{spread.level.plot} is replaced by \code{\link{spreadLevelPlot}}.
}
