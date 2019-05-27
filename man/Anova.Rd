%-------------------------------------------------------------------------------------
% Revision history:
% checked in 2008-12-29 by J. Fox (corresponds to version 1.2-10 of car)
%	2009-01-16 updated doc to correspond to changes in linearHypothesis. J. Fox
%   2009-09-16 updated to reflect new singular.ok argument
%   2009-09-26 updated for removal from car. J. Fox
%   2009-12-22 updated to reflect new imatrix argument to Anova.mlm(). J. Fox
%   2012-02-28 updated to reflect new test.statistic argument to Anova.mer(). J. Fox
%   2012-03-01 removed commented-out output listings. J. Fox
%   2013-06-17 update for summary.Anova.mlm and print.summary.Anova.mlm. J. Fox
%   2013-06-20 added Anova.merMod(). J. Fox
%   2014-08-09: added vcov. argument to Anova.lm(). J. Fox
%   2014-09-23: added Anova.rlm(). J. Fox
%   2015-02-18: updated info about when Anova.default() works. J. Fox
%   2015-09-04: added Anova.coxme(). J. Fox
%   2016-06-03: added SSP and SSPE arguments to print.summary.ANova.mlm(). J. Fox
%   2016-06-25: provision to print univariate ANOVAs for a mlm. J. Fox
%   2017-11-02: added warning about KR F-tests
%   2017-11-07: make compatible with vcov() in re-devel. J. Fox
%-------------------------------------------------------------------------------------

\name{Anova}
\alias{Anova}
\alias{Anova.lm}
\alias{Anova.lme}
\alias{Anova.aov}
\alias{Anova.glm}
\alias{Anova.multinom}
\alias{Anova.polr}
\alias{Anova.mer}
\alias{Anova.merMod}
\alias{Anova.mlm}
\alias{Anova.manova}
\alias{Manova}
\alias{Manova.mlm}
\alias{print.Anova.mlm}
\alias{summary.Anova.mlm}
\alias{print.summary.Anova.mlm}
\alias{print.univaov}
\alias{as.data.frame.univaov}
\alias{Anova.coxph}
\alias{Anova.svyglm}
\alias{Anova.rlm}
\alias{Anova.coxme}
\alias{Anova.default}
\title{Anova Tables for Various Statistical Models}
\description{
  Calculates type-II or type-III analysis-of-variance tables for
  model objects produced by \code{lm}, \code{glm}, \code{multinom}
  (in the \pkg{nnet} package), \code{polr} (in the \pkg{MASS}
  package), \code{coxph} (in the \pkg{survival} package),
  \code{coxme} (in the \pkg{coxme} pckage),
  \code{svyglm} (in the \pkg{survey} package), \code{rlm} (in the \pkg{MASS} package),
  \code{lmer} in the \pkg{lme4} package,
  \code{lme} in the \pkg{nlme} package, and (by the default method) for most
  models with a linear predictor and asymptotically normal coefficients (see details below). For linear
  models, F-tests are calculated; for generalized linear models,
  likelihood-ratio chisquare, Wald chisquare, or F-tests are calculated;
  for multinomial logit and proportional-odds logit models, likelihood-ratio
  tests are calculated.  Various test statistics are provided for multivariate
  linear models produced by \code{lm} or \code{manova}. Partial-likelihood-ratio tests
  or Wald tests are provided for Cox models. Wald chi-square tests are provided for fixed effects in
  linear and generalized linear mixed-effects models. Wald chi-square or F tests are provided
  in the default case.
}
\usage{
Anova(mod, ...)

Manova(mod, ...)

\method{Anova}{lm}(mod, error, type=c("II","III", 2, 3),
	white.adjust=c(FALSE, TRUE, "hc3", "hc0", "hc1", "hc2", "hc4"),
	vcov.=NULL, singular.ok, ...)

\method{Anova}{aov}(mod, ...)

\method{Anova}{glm}(mod, type=c("II","III", 2, 3),
    test.statistic=c("LR", "Wald", "F"),
    error, error.estimate=c("pearson", "dispersion", "deviance"),
    singular.ok, ...)

\method{Anova}{multinom}(mod, type = c("II","III", 2, 3), ...)

\method{Anova}{polr}(mod, type = c("II","III", 2, 3), ...)

\method{Anova}{mlm}(mod, type=c("II","III", 2, 3), SSPE, error.df,
    idata, idesign, icontrasts=c("contr.sum", "contr.poly"), imatrix,
    test.statistic=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),...)

\method{Anova}{manova}(mod, ...)

\method{Manova}{mlm}(mod, ...)

\method{print}{Anova.mlm}(x, ...)

\method{summary}{Anova.mlm}(object, test.statistic, univariate=object$repeated,
    multivariate=TRUE, p.adjust.method, ...)

\method{print}{summary.Anova.mlm}(x, digits = getOption("digits"),
    SSP=TRUE, SSPE=SSP, ... )

\method{print}{univaov}(x, digits = max(getOption("digits") - 2L, 3L),
                          style=c("wide", "long"),
                          by=c("response", "term"),
                          ...)

\method{as.data.frame}{univaov}(x, row.names, optional, by=c("response", "term"), ...)

\method{Anova}{coxph}(mod, type=c("II", "III", 2, 3),
	test.statistic=c("LR", "Wald"), ...)

\method{Anova}{coxme}(mod, type=c("II", "III", 2, 3),
    test.statistic=c("Wald", "LR"), ...)

\method{Anova}{lme}(mod, type=c("II","III", 2, 3),
		vcov.=vcov(mod, complete=FALSE), singular.ok, ...)

\method{Anova}{mer}(mod, type=c("II", "III", 2, 3),
	test.statistic=c("Chisq", "F"), vcov.=vcov(mod, complete=FALSE), singular.ok, ...)

\method{Anova}{merMod}(mod, type=c("II", "III", 2, 3),
    test.statistic=c("Chisq", "F"), vcov.=vcov(mod, complete=FALSE), singular.ok, ...)

\method{Anova}{svyglm}(mod, ...)

\method{Anova}{rlm}(mod, ...)

\method{Anova}{default}(mod, type=c("II", "III", 2, 3),
	test.statistic=c("Chisq", "F"), vcov.=vcov(mod, complete=FALSE),
	singular.ok, ...)
}

\arguments{
  \item{mod}{\code{lm}, \code{aov}, \code{glm}, \code{multinom}, \code{polr}
    \code{mlm}, \code{coxph}, \code{coxme}, \code{lme}, \code{mer}, \code{merMod}, \code{svyglm},
    \code{rlm}, or other suitable model object.}
  \item{error}{for a linear model, an \code{lm} model object from which the
    error sum of squares and degrees of freedom are to be calculated. For
    F-tests for a generalized linear model, a \code{glm} object from which the
    dispersion is to be estimated. If not specified, \code{mod} is used.}
  \item{type}{type of test, \code{"II"}, \code{"III"}, \code{2}, or \code{3}.  Roman numerals are equivalent to
    the corresponding Arabic numerals.}
  \item{singular.ok}{defaults to \code{TRUE} for type-II tests, and \code{FALSE}
    for type-III tests where the tests for models with aliased coefficients
    will not be straightforwardly interpretable;
    if \code{FALSE}, a model with aliased coefficients produces an error.}
  \item{test.statistic}{for a generalized linear model, whether to calculate
    \code{"LR"} (likelihood-ratio), \code{"Wald"}, or \code{"F"} tests; for a Cox
    or Cox mixed-effects model, whether to calculate \code{"LR"} (partial-likelihood ratio) or
    \code{"Wald"} tests; in the default case or for linear mixed models fit by
    \code{lmer}, whether to calculate Wald \code{"Chisq"} or Kenward-Roger
    \code{"F"} tests with Satterthwaite degrees of freedom (\emph{warning:} the KR F-tests
    can be very time-consuming).
    For a multivariate linear model, the multivariate test statistic to compute --- one of
    \code{"Pillai"}, \code{"Wilks"}, \code{"Hotelling-Lawley"}, or \code{"Roy"},
    with \code{"Pillai"} as the default. The \code{summary} method for \code{Anova.mlm}
    objects permits the specification of more than one multivariate
    test statistic, and the default is to report all four.}
  \item{error.estimate}{for F-tests for a generalized linear model, base the
    dispersion estimate on the Pearson residuals (\code{"pearson"}, the default); use the
    dispersion estimate in the model object (\code{"dispersion"}); or base the dispersion estimate on
    the residual deviance (\code{"deviance"}). For binomial or Poisson GLMs, where the dispersion
    is fixed to 1, setting \code{error.estimate="dispersion"} is changed to \code{"pearson"},
    with a warning.}
  \item{white.adjust}{if not \code{FALSE}, the default,
  	tests use a heteroscedasticity-corrected coefficient
    covariance matrix; the various values of the argument specify different corrections.
    See the documentation for \code{\link{hccm}} for details. If \code{white.adjust=TRUE}
    then the \code{"hc3"} correction is selected.}
  \item{SSPE}{For \code{Anova} for a multivariate linear model, the
    error sum-of-squares-and-products matrix; if missing, will be computed
    from the residuals of the model; for the \code{print} method for the \code{summary} of
    an \code{Anova} of a multivariate linear model,
    whether or not to print the error SSP matrix (defaults to \code{TRUE}).}
  \item{SSP}{if \code{TRUE} (the default), print the sum-of-squares and
    cross-products matrix for the hypothesis and the response-transformation matrix.}
  \item{error.df}{The degrees of freedom for error; if missing, will be taken from the model.}
  \item{idata}{an optional data frame giving a factor or factors defining the
    intra-subject model for multivariate repeated-measures data. See
    \emph{Details} for an explanation of the intra-subject design and for
    further explanation of the other arguments relating to intra-subject factors.}
  \item{idesign}{a one-sided model formula using the ``data'' in \code{idata} and
    specifying the intra-subject design.}
  \item{icontrasts}{names of contrast-generating functions to be applied by default
    to factors and ordered factors, respectively, in the within-subject
    ``data''; the contrasts must produce an intra-subject model
    matrix in which different terms are orthogonal. The default is
    \code{c("contr.sum", "contr.poly")}.}
  \item{imatrix}{as an alternative to specifying \code{idata}, \code{idesign}, and
  	(optionally) \code{icontrasts}, the model matrix for the within-subject design
  	can be given directly in the form of list of named elements. Each element gives
  	the columns of the within-subject model matrix for a term to be tested, and must
  	have as many rows as there are responses; the columns of the within-subject model
  	matrix for different terms must be mutually orthogonal.}
  \item{x, object}{object of class \code{"Anova.mlm"} to print or summarize.}
  \item{multivariate, univariate}{compute and print multivariate and univariate tests for a repeated-measures
    ANOVA or multivariate linear model; the default is \code{TRUE} for both for repeated measures and \code{TRUE}
    for \code{multivariate} for a multivariate linear model.}
  \item{p.adjust.method}{if given for a multivariate linear model when univariate tests are requested, the
    univariate tests are corrected for simultaneous inference by term; if specified, should be one of the methods
    recognized by \code{\link{p.adjust}} or \code{TRUE}, in which case the default (Holm) adjustment is used.}
  \item{digits}{minimum number of significant digits to print.}
  \item{style}{for printing univariate tests if requested for a multivariate linear model; one of \code{"wide"},
    the default, or \code{"long"}.}
  \item{by}{if univariate tests are printed in \code{"long"} \code{style}, they can be ordered \code{by}
    \code{"response"}, the default, or by \code{"term"}.}
  \item{row.names, optional}{not used.}
  \item{vcov.}{in the \code{default} method, an optional coefficient-covariance matrix or function
    to compute a covariance matrix, computed by default by applying the generic \code{vcov} function to the model object.
    A similar argument may be supplied to the \code{lm} method, and the default (\code{NULL}) is to ignore the argument;
    if both \code{vcov.} and \code{white.adjust} are supplied to the \code{lm} method, the latter is used.}
  \item{\dots}{do not use.}
}
\details{
  The designations "type-II" and "type-III" are borrowed from SAS, but the
  definitions used here do not correspond precisely to those employed by SAS.
  Type-II tests are calculated according to the principle of marginality,
  testing each term after all others, except ignoring the term's higher-order relatives;
  so-called type-III tests violate marginality, testing
  each term in the model after all of the others. This definition of Type-II tests
  corresponds to the tests produced by SAS for analysis-of-variance models, where all of the predictors
  are factors, but not more generally (i.e., when there are quantitative predictors).
  Be very careful in formulating the model for type-III tests, or the hypotheses tested
  will not make sense.

  As implemented here, type-II Wald tests are a generalization of the linear hypotheses used to generate
  these tests in linear models.

  For tests for linear models, multivariate linear models, and Wald tests for generalized linear models,
  Cox models, mixed-effects models, generalized linear models fit to survey data, and in the default case,
  \code{Anova} finds the test statistics without refitting the model. The \code{svyglm} method simply
  calls the \code{default} method and therefore can take the same arguments.

  The standard R \code{anova} function calculates sequential ("type-I") tests.
  These rarely test interesting hypotheses in unbalanced designs.

  A MANOVA for a multivariate linear model (i.e., an object of
  class \code{"mlm"} or \code{"manova"}) can optionally include an
  intra-subject repeated-measures design.
  If the intra-subject design is absent (the default), the multivariate
  tests concern all of  the response variables.
  To specify a repeated-measures design, a data frame is provided defining the repeated-measures factor or
  factors
  via \code{idata}, with default contrasts given by the \code{icontrasts}
  argument. An intra-subject model-matrix is generated from the formula
  specified by the \code{idesign} argument; columns of the model matrix
  corresponding to different terms in the intra-subject model must be orthogonal
  (as is insured by the default contrasts). Note that the contrasts given in
  \code{icontrasts} can be overridden by assigning specific contrasts to the
  factors in \code{idata}. As an alternative, the within-subjects model matrix
  can be specified directly via the \code{imatrix} argument.
  \code{Manova} is essentially a synonym for \code{Anova}
  for multivariate linear models.

  If univariate tests are requested for the \code{summary} of a multivariate linear model, the object returned
  contains a \code{univaov} component of \code{"univaov"}; \code{print} and \code{as.data.frame} methods are
  provided for the \code{"univaov"} class.

  For the default method to work, the model object must contain a standard
  \code{terms} element, and must respond to the \code{vcov}, \code{coef}, and \code{model.matrix} functions.
  If any of these requirements is missing, then it may be possible to supply it reasonably simply (e.g., by
  writing a missing \code{vcov} method for the class of the model object).
}
\value{
  An object of class \code{"anova"}, or \code{"Anova.mlm"}, which usually is printed.
  For objects of class \code{"Anova.mlm"}, there is also a \code{summary} method,
  which provides much more detail than the \code{print} method about the MANOVA, including
  traditional mixed-model univariate F-tests with Greenhouse-Geisser and Huynh-Feldt
  corrections.
}
\references{
  Fox, J. (2016)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Third Edition. Sage.

  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.

  Hand, D. J., and Taylor, C. C. (1987)
  \emph{Multivariate Analysis of Variance and Repeated Measures: A Practical
  Approach for Behavioural Scientists.} Chapman and Hall.

  O'Brien, R. G., and Kaiser, M. K. (1985)
  MANOVA method for analyzing repeated measures designs: An extensive primer.
  \emph{Psychological Bulletin} \bold{97}, 316--333.
}
\author{John Fox \email{jfox@mcmaster.ca}; the code for the Mauchly test and Greenhouse-Geisser and Huynh-Feldt
corrections for non-spericity in repeated-measures ANOVA are adapted from the functions
\code{stats:::stats:::mauchly.test.SSD} and \code{stats:::sphericity} by R Core; \code{summary.Anova.mlm}
and \code{print.summary.Anova.mlm} incorporates code contributed by Gabriel Baud-Bovy.}

\section{Warning}{
Be careful of type-III tests: For a traditional multifactor ANOVA model with interactions, for example, these tests will normally only be sensible when using contrasts that, for different terms, are
orthogonal in the row-basis of the model, such as those produced by \code{\link{contr.sum}}, \code{\link{contr.poly}}, or \code{\link{contr.helmert}}, but \emph{not} by the default
\code{\link{contr.treatment}}. In a model that contains factors, numeric covariates, and interactions, main-effect tests for factors will be for differences over the origin. In contrast (pun intended),
type-II tests are invariant with respect to (full-rank) contrast coding. If you don't understand this issue, then you probably shouldn't use \code{Anova} for type-III tests.
}

\seealso{ \code{\link{linearHypothesis}}, \code{\link[stats]{anova}}
\code{\link[stats]{anova.lm}}, \code{\link[stats]{anova.glm}},
\code{\link[stats]{anova.mlm}}, \code{\link[survival]{anova.coxph}}, \code{\link[survey]{svyglm}}.}

\examples{

## Two-Way Anova

mod <- lm(conformity ~ fcategory*partner.status, data=Moore,
  contrasts=list(fcategory=contr.sum, partner.status=contr.sum))
Anova(mod)
Anova(mod, type=3)  # note use of contr.sum in call to lm()

## One-Way MANOVA
## See ?Pottery for a description of the data set used in this example.

summary(Anova(lm(cbind(Al, Fe, Mg, Ca, Na) ~ Site, data=Pottery)))

## MANOVA for a randomized block design (example courtesy of Michael Friendly:
##  See ?Soils for description of the data set)

soils.mod <- lm(cbind(pH,N,Dens,P,Ca,Mg,K,Na,Conduc) ~ Block + Contour*Depth,
    data=Soils)
Manova(soils.mod)
summary(Anova(soils.mod), univariate=TRUE, multivariate=FALSE,
    p.adjust.method=TRUE)

## a multivariate linear model for repeated-measures data
## See ?OBrienKaiser for a description of the data set used in this example.

phase <- factor(rep(c("pretest", "posttest", "followup"), c(5, 5, 5)),
    levels=c("pretest", "posttest", "followup"))
hour <- ordered(rep(1:5, 3))
idata <- data.frame(phase, hour)
idata

mod.ok <- lm(cbind(pre.1, pre.2, pre.3, pre.4, pre.5,
                     post.1, post.2, post.3, post.4, post.5,
                     fup.1, fup.2, fup.3, fup.4, fup.5) ~  treatment*gender,
                data=OBrienKaiser)
(av.ok <- Anova(mod.ok, idata=idata, idesign=~phase*hour))

summary(av.ok, multivariate=FALSE)

## A "doubly multivariate" design with two  distinct repeated-measures variables
## (example courtesy of Michael Friendly)
## See ?WeightLoss for a description of the dataset.

imatrix <- matrix(c(
	1,0,-1, 1, 0, 0,
	1,0, 0,-2, 0, 0,
	1,0, 1, 1, 0, 0,
	0,1, 0, 0,-1, 1,
	0,1, 0, 0, 0,-2,
	0,1, 0, 0, 1, 1), 6, 6, byrow=TRUE)
colnames(imatrix) <- c("WL", "SE", "WL.L", "WL.Q", "SE.L", "SE.Q")
rownames(imatrix) <- colnames(WeightLoss)[-1]
(imatrix <- list(measure=imatrix[,1:2], month=imatrix[,3:6]))
contrasts(WeightLoss$group) <- matrix(c(-2,1,1, 0,-1,1), ncol=2)
(wl.mod<-lm(cbind(wl1, wl2, wl3, se1, se2, se3)~group, data=WeightLoss))
Anova(wl.mod, imatrix=imatrix, test="Roy")

## mixed-effects models examples:

\dontrun{
	library(nlme)
	example(lme)
	Anova(fm2)
}

\dontrun{
	library(lme4)
	example(glmer)
	Anova(gm1)
}

}
\keyword{htest}
\keyword{models}
\keyword{regression}
