\name{effect}
\alias{effect}
\alias{all.effects}
\alias{print.effect}
\alias{print.effect.list}
\alias{summary.effect}
\alias{summary.effect.list}
\alias{plot.effect}
\alias{plot.effect.list}
\alias{as.data.frame.effect}

\title{Functions For Constructing Effect Plots}
\description{
  \code{effect} constructs an \code{"effect"} object for a high-order term in a
  linear or generalized linear model, absorbing the lower-order terms marginal
  to the term in question, and averaging over other terms in the model. 
  
  \code{all.effects} identifies all of the high-order terms in a model and returns
  a list of \code{"effect"} objects (i.e., an object of type \code{"effect.list"}).
}

\usage{
effect(term, mod, xlevels=list(), default.levels=10, se=TRUE, 
    confidence.level=.95, transformation=family(mod)$linkinv, typical=mean)
    
all.effects(mod, ...)

\method{print}{effect}(x, ...)
%% \method{print}{effect.list}(x, ...)
\method{summary}{effect}(object, ...)
%% \method{summary{effect.list}(object, ...)
\method{plot}{effect}(x, x.var=which.max(levels), 
    z.var=which.min(levels), multiline=is.null(x$se), rug=TRUE,
    ylab=x$response, colors=palette(), symbols=1:10, lines=1:10,
    factor.names=TRUE, ...)
\method{plot}{effect.list}(x, selection, ...)
as.data.frame.effect(x, row.names=NULL, optional=TRUE)
}

\arguments{
  \item{term}{the quoted name of a high-order term in the model.}
  \item{mod}{an object of class \code{"lm"} or \code{"glm"}.}
  \item{xlevels}{an optional list of values at which to set covariates,
    with components of the form \code{covariate.name = vector.of.values}.}
  \item{default.levels}{number of values for covariates that are not
    specified explicitly via \code{xlevels}; covariate values set by
    default are evenly spaced between the minimum and maximum values in
    the data.}
  \item{se}{if \code{TRUE}, the default, calculate standard errors and
    confidence limits for the effects.}
  \item{confidence.level}{level at which to compute confidence limits
    based on the standard-normal distribution; the default is \code{0.95}.}
  \item{transformation}{a function to be applied to the calculated effects
    and confidence limits; the default is the inverse-link function for the
    model, which is a the identity function for a linear model
    (i.e., no transformation.)}
  \item{typical}{a function to be applied to the columns of the model matrix
    over which the effect is "averaged"; the default is \code{mean}.}
  \item{x, object}{an object of type \code{"effect"} or \code{"effect.list"},
    as appropriate.}
  \item{x.var}{the index (number) of the covariate or factor to place on the 
    horizontal axis of each panel of the effect plot. The default is the 
    predictor with the largest number of levels or values.}
  \item{multiline}{if \code{TRUE}, each panel of the display represents combinations
    of values of two predictors, with one predictor (corresponding to \code{x.var})
    on the horzontal axis, and the other (corresponding to \code{z.var}) used to define
    lines in the graph; defaults to \code{TRUE} if there are no standard errors in
    the object being plotted, and \code{FALSE} otherwise.}
  \item{rug}{if \code{TRUE}, the default, a rug plot is shown giving the marginal 
    distribution of the predictor on the horizontal axis, if this predictor is
    a covariate.}
  \item{z.var}{the index (number) of the covariate or factor for which
    individual lines are to be drawn in each panel of the effect plot. The default is the 
    predictor with the smallest number of levels or values. This argument is only
    used if \code{multiline = TRUE}.}
  \item{colors}{\code{colors[1]} is used to plot effects, \code{colors[2]} to plot 
    confidence bands. In a mulitline plot, the successive \code{colors} correspond
    to the levels of the \code{z.var} covariate or factor.}
  \item{symbols, lines}{corresponding to the levels of the \code{z.var} covariate
    or factor on a multiline plot. These arguments are used only if \code{multiline = TRUE};
    in this case a legend is drawn at the top of the display.}
  \item{factor.names}{a logical value, default \code{TRUE}, that controls the inclusion of
    factor names in conditioning-variable labels.}
  \item{ylab}{the label for the vertical axis of the effect plot; the default is the
    response variable for the model from which the effect was computed.}
  \item{selection}{the optional index (number) or quoted name of the effect in an effect
    list to be plotted; if not supplied, a menu of high-order terms is presented.}
  \item{...}{arguments to be passed down.}
  \item{row.names, optional}{not used.}
}

\details{
  Normally, the functions to be used directly are \code{all.effects}, to return
  a list of effects, and the generic \code{plot} function to plot the effects.
  Plots are drawn using the \code{\link[lattice]{xyplot}} function in the 
  \code{lattice} package. Effects may also be printed (implicitly or explicitly via
  \code{print}) or summarized (using \code{summary}).
  
  In calculating effects, the strategy for "safe" prediction described
  in Hastie (1992: Sec. 7.3.3) is employed.
}

\value{
  \code{effect} returns a list with the following components:
  \item{term}{the high-order term to which the effect pertains.}
  \item{formula}{the complete model formula.}
  \item{response}{a character string giving the response variable.}
  \item{variables}{a list with information about each predictor, including
    its name, whether it is a factor, and its levels or values.}
  \item{effect}{a one-column matrix of fitted values, representing the effect
    on the original scale; this is a ravelled table, representing all
    combinations of predictor values.}
  \item{fit}{a one-column matrix of fitted values, representing the effect on
    the transformed scale.}
  \item{x}{a data frame, the columns of which are the predictors, and the
    rows of which give all combinations of values of the predictors.}
  \item{model.matrix}{the model matrix from which the effect was calculated.}
  \item{data} {a data frame with the data on which the fitted model was based.}
  \item{discrepancy}{the percentage discrepancy for the 'safe' predictions of
    the original fit; should be very close to 0.}
  \item{se}{a vector or standard errors for the effect, on the original scale.}
  \item{lower, upper}{one-column matrices of confidence limits, on the
    transformed scale.}
  \item{confidence.level}{corresponding to the confidence limits.}
}

\references{
    Fox, J. (1987)
    Effect displays for generalized linear models.
    \emph{Sociological Methodology}
    \bold{17}, 347--361.
    
    Hastie, T. J. (1992)
    Generalized additive models.
    In Chambers, J. M., and Hastie, T. J. (eds.)
    \emph{Statistical Models in S}, Wadsworth.
    }
    
    
    
\author{John Fox jfox@mcmaster.ca}


\seealso{\code{\link[lattice]{xyplot}}}

\examples{
data(Cowles)
mod.cowles <- glm(volunteer ~ sex + neuroticism*extraversion, 
    data=Cowles, family=binomial)
eff.cowles <- all.effects(mod.cowles, xlevels=list(neuroticism=0:24, 
    extraversion=seq(0, 24, 6)))
eff.cowles
    \dontrun{
    model: volunteer ~ sex + neuroticism * extraversion
    
    sex effect
    sex
    female      male 
    0.4409441 0.3811941 
    
    neuroticism*extraversion effect
            extraversion
    neuroticism          0         6        12        18        24
            0  0.07801066 0.1871263 0.3851143 0.6301824 0.8225756
            1  0.08636001 0.1963396 0.3870453 0.6200668 0.8083638
            2  0.09551039 0.2058918 0.3889798 0.6098458 0.7932997
            3  0.10551835 0.2157839 0.3909179 0.5995275 0.7773775
           . . .
            23 0.51953129 0.4747277 0.4303273 0.3870199 0.3454282
            24 0.54709527 0.4895731 0.4323256 0.3768303 0.3243880
    }
plot(eff.cowles, 'sex', ylab="Prob(Volunteer)")
    \dontrun{
    Loading required package: lattice 
    Loading required package: grid 
    }
plot(eff.cowles, 'neuroticism:extraversion', ylab="Prob(Volunteer)")
plot(eff.cowles, 'neuroticism:extraversion', multiline=TRUE, 
    ylab="Prob(Volunteer)")
    
plot(effect('sex:neuroticism:extraversion',
    glm(volunteer ~ sex*neuroticism*extraversion, data=Cowles, family=binomial),
    xlevels=list(neuroticism=0:24, extraversion=seq(0, 24, 6))), multiline=TRUE)
    
\dontrun{
data(Prestige)
mod.pres <- lm(prestige ~ log(income, 10) + poly(education, 3) + poly(women, 2), 
    data=Prestige)
eff.pres <- all.effects(mod.pres, default.levels=50)
plot(eff.pres)
}

}

\keyword{hplot}
\keyword{models}
