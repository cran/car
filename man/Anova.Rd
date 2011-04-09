%-------------------------------------------------------------------------------
% Revision history:
% checked in 2008-12-29 by J. Fox (corresponds to version 1.2-10 of car)
%	2009-01-16 updated doc to correspond to changes in linearHypothesis. J. Fox
%   2009-09-16 updated to reflect new singular.ok argument
%   2009-09-26 updated for removal from car. J. Fox
%   2009-12-22 updated to reflect new imatrix argument to Anova.mlm(). J. Fox
%-------------------------------------------------------------------------------

\name{Anova}
\alias{Anova}
\alias{Anova.lm}
\alias{Anova.lme}
\alias{Anova.aov}
\alias{Anova.glm}
\alias{Anova.multinom}
\alias{Anova.polr}
\alias{Anova.mer}
\alias{Anova.mlm}
\alias{Anova.manova}
\alias{Manova}
\alias{Manova.mlm}
\alias{print.Anova.mlm}
\alias{summary.Anova.mlm}
\alias{Anova.coxph}
\alias{Anova.default}
\title{Anova Tables for Various Statistical Models}
\description{
  Calculates type-II or type-III analysis-of-variance tables for
  model objects produced by \code{lm}, \code{glm}, \code{multinom} 
  (in the \pkg{nnet} package), \code{polr} (in the \pkg{MASS}
  package), \code{coxph} (in the \pkg{survival} package), \code{lmer} in the \pkg{lme4} package,
  \code{lme} in the \pkg{nlme} package, and for any
  model with a linear predictor and asymptotically normal coefficients that
  responds to the \code{vcov} and \code{coef} functions.  For linear
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
	singular.ok, ...)

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

\method{summary}{Anova.mlm}(object, test.statistic, multivariate=TRUE, 
    univariate=TRUE, digits=getOption("digits"), ...)
    
\method{Anova}{coxph}(mod, type=c("II","III", 2, 3), 
	test.statistic=c("LR", "Wald"), ...)
	
\method{Anova}{lme}(mod, type=c("II","III", 2, 3),
		vcov.=vcov(mod), singular.ok, ...)

\method{Anova}{mer}(mod, type=c("II","III", 2, 3),
		vcov.=vcov(mod), singular.ok, ...)
	
\method{Anova}{default}(mod, type=c("II","III", 2, 3), 
	test.statistic=c("Chisq", "F"), vcov.=vcov(mod), 
	singular.ok, ...)
}

\arguments{
  \item{mod}{\code{lm}, \code{aov}, \code{glm}, \code{multinom}, \code{polr}
    \code{mlm}, \code{coxph} or other suitable model object.}
  \item{error}{for a linear model, an \code{lm} model object from which the
    error sum of squares and degrees of freedom are to be calculated. For 
    F-tests for a generalized linear model, a \code{glm} object from which the
    dispersion is to be estimated. If not specified, \code{mod} is used.}
  \item{type}{type of test, \code{"II"}, \code{"III"}, \code{2}, or \code{3}.}
  \item{singular.ok}{defaults to \code{TRUE} for type-II tests, and \code{FALSE}
    for type-III tests (where the tests for models with aliased coefficients
    will not be straightforwardly interpretable); 
    if \code{FALSE}, a model with aliased coefficients produces an error.}
  \item{test.statistic}{for a generalized linear model, whether to calculate 
    \code{"LR"} (likelihood-ratio), \code{"Wald"}, or \code{"F"} tests; for a Cox
    model, whether to calculate \code{"LR"} (partial-likelihood ratio) or
    \code{"Wald"} tests; in the default case, whether to calculate Wald \code{"Chisq"} or
    \code{"F"} tests.
    For a multivariate linear model, the multivariate test statistic to compute --- one of
    \code{"Pillai"}, \code{"Wilks"}, \code{"Hotelling-Lawley"}, or \code{"Roy"}, 
    with \code{"Pillai"} as the default. The \code{summary} method for \code{Anova.mlm}
    objects permits the specification of more than one multivariate
    test statistic, and the default is to report all four.}
  \item{error.estimate}{for F-tests for a generalized linear model, base the
    dispersion estimate on the Pearson residuals (\code{"pearson"}, the default); use the
    dispersion estimate in the model object (\code{"dispersion"}), which, e.g., is
    fixed to 1 for binomial and Poisson models; or base the dispersion estimate on
    the residual deviance (\code{"deviance"}).}
  \item{white.adjust}{if not \code{FALSE}, the default, 
  	tests use a heteroscedasticity-corrected coefficient
    covariance matrix; the various values of the argument specify different corrections.
    See the documentation for \code{\link{hccm}} for details. If \code{white.adjust=TRUE}
    then the \code{"hc3"} correction is selected.} 
  \item{SSPE}{The error sum-of-squares-and-products matrix; if missing, will be computed
    from the residuals of the model.}
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
  \item{multivariate, univariate}{print multivariate and univariate tests for a repeated-measures
    ANOVA; the default is \code{TRUE} for both.}
  \item{digits}{minimum number of significant digits to print.}
  \item{vcov.}{an optional coefficient-covariance matrix, computed by default by applying the 
  	generic \code{vcov} function to the model object.}
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
  Cox models, and in the default case,
  \code{Anova} finds the test statistics without refitting the model.
  
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
}
\value{
  An object of class \code{"anova"}, or \code{"Anova.mlm"}, which usually is printed.
  For objects of class \code{"Anova.mlm"}, there is also a \code{summary} method, 
  which provides much more detail than the \code{print} method about the MANOVA, including
  traditional mixed-model univariate F-tests with Greenhouse-Geisser and Huynh-Feldt
  corrections.
}
\references{ 
  Fox, J. (2008)
  \emph{Applied Regression Analysis and Generalized Linear Models},
  Second Edition. Sage.  
  
  Fox, J. and Weisberg, S. (2011) 
  \emph{An R Companion to Applied Regression}, Second Edition, Sage.
  
  Hand, D. J., and Taylor, C. C. (1987)
  \emph{Multivariate Analysis of Variance and Repeated Measures: A Practical
  Approach for Behavioural Scientists.} Chapman and Hall.
  
  O'Brien, R. G., and Kaiser, M. K. (1985)
  MANOVA method for analyzing repeated measures designs: An extensive primer.
  \emph{Psychological Bulletin} \bold{97}, 316--333.
}
\author{John Fox \email{jfox@mcmaster.ca}}

\section{Warning}{Be careful of type-III tests.}

\seealso{ \code{\link{linearHypothesis}}, \code{\link[stats]{anova}}
\code{\link[stats]{anova.lm}}, \code{\link[stats]{anova.glm}}, 
\code{\link[stats]{anova.mlm}}, \code{\link[survival]{anova.coxph}}.}

\examples{

## Two-Way Anova

mod <- lm(conformity ~ fcategory*partner.status, data=Moore, 
  contrasts=list(fcategory=contr.sum, partner.status=contr.sum))
Anova(mod)
## Anova Table (Type II tests)
##
## Response: conformity
##                         Sum Sq Df F value   Pr(>F)
## fcategory                 11.61  2  0.2770 0.759564
## partner.status           212.21  1 10.1207 0.002874
## fcategory:partner.status 175.49  2  4.1846 0.022572
## Residuals                817.76 39                 
Anova(mod, type="III")
## Anova Table (Type III tests)
##
## Response: conformity
##                          Sum Sq Df  F value    Pr(>F)
## (Intercept)              5752.8  1 274.3592 < 2.2e-16
## fcategory                  36.0  2   0.8589  0.431492
## partner.status            239.6  1  11.4250  0.001657
## fcategory:partner.status  175.5  2   4.1846  0.022572
## Residuals                 817.8 39

## One-Way MANOVA
## See ?Pottery for a description of the data set used in this example.

summary(Anova(lm(cbind(Al, Fe, Mg, Ca, Na) ~ Site, data=Pottery)))

## Type II MANOVA Tests:
## 
## Sum of squares and products for error:
##            Al          Fe          Mg          Ca         Na
## Al 48.2881429  7.08007143  0.60801429  0.10647143 0.58895714
## Fe  7.0800714 10.95084571  0.52705714 -0.15519429 0.06675857
## Mg  0.6080143  0.52705714 15.42961143  0.43537714 0.02761571
## Ca  0.1064714 -0.15519429  0.43537714  0.05148571 0.01007857
## Na  0.5889571  0.06675857  0.02761571  0.01007857 0.19929286
## 
## ------------------------------------------
##  
## Term: Site 
## 
## Sum of squares and products for the hypothesis:
##             Al          Fe          Mg         Ca         Na
## Al  175.610319 -149.295533 -130.809707 -5.8891637 -5.3722648
## Fe -149.295533  134.221616  117.745035  4.8217866  5.3259491
## Mg -130.809707  117.745035  103.350527  4.2091613  4.7105458
## Ca   -5.889164    4.821787    4.209161  0.2047027  0.1547830
## Na   -5.372265    5.325949    4.710546  0.1547830  0.2582456
## 
## Multivariate Tests: Site
##                        Df test stat  approx F   num Df   den Df     Pr(>F)    
## Pillai            3.00000   1.55394   4.29839 15.00000 60.00000 2.4129e-05 ***
## Wilks             3.00000   0.01230  13.08854 15.00000 50.09147 1.8404e-12 ***
## Hotelling-Lawley  3.00000  35.43875  39.37639 15.00000 50.00000 < 2.22e-16 ***
## Roy               3.00000  34.16111 136.64446  5.00000 20.00000 9.4435e-15 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

## MANOVA for a randomized block design (example courtesy of Michael Friendly:
##  See ?Soils for description of the data set)

soils.mod <- lm(cbind(pH,N,Dens,P,Ca,Mg,K,Na,Conduc) ~ Block + Contour*Depth, 
    data=Soils)
Manova(soils.mod)

## Type II MANOVA Tests: Pillai test statistic
##                Df test stat approx F num Df den Df    Pr(>F)    
## Block           3    1.6758   3.7965     27     81 1.777e-06 ***
## Contour         2    1.3386   5.8468     18     52 2.730e-07 ***
## Depth           3    1.7951   4.4697     27     81 8.777e-08 ***
## Contour:Depth   6    1.2351   0.8640     54    180    0.7311    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 


## a multivariate linear model for repeated-measures data
## See ?OBrienKaiser for a description of the data set used in this example.

phase <- factor(rep(c("pretest", "posttest", "followup"), c(5, 5, 5)),
    levels=c("pretest", "posttest", "followup"))
hour <- ordered(rep(1:5, 3))
idata <- data.frame(phase, hour)
idata
##       phase hour
## 1   pretest    1
## 2   pretest    2
## 3   pretest    3
## 4   pretest    4
## 5   pretest    5
## 6  posttest    1
## 7  posttest    2
## 8  posttest    3
## 9  posttest    4
## 10 posttest    5
## 11 followup    1
## 12 followup    2
## 13 followup    3
## 14 followup    4
## 15 followup    5

mod.ok <- lm(cbind(pre.1, pre.2, pre.3, pre.4, pre.5, 
                     post.1, post.2, post.3, post.4, post.5, 
                     fup.1, fup.2, fup.3, fup.4, fup.5) ~  treatment*gender, 
                data=OBrienKaiser)
(av.ok <- Anova(mod.ok, idata=idata, idesign=~phase*hour)) 
## Type II Repeated Measures MANOVA Tests: Pillai test statistic
##                             Df test stat approx F num Df den Df    Pr(>F)    
## treatment                    2    0.4809   4.6323      2     10 0.0376868 *  
## gender                       1    0.2036   2.5558      1     10 0.1409735    
## treatment:gender             2    0.3635   2.8555      2     10 0.1044692    
## phase                        1    0.8505  25.6053      2      9 0.0001930 ***
## treatment:phase              2    0.6852   2.6056      4     20 0.0667354 .  
## gender:phase                 1    0.0431   0.2029      2      9 0.8199968    
## treatment:gender:phase       2    0.3106   0.9193      4     20 0.4721498    
## hour                         1    0.9347  25.0401      4      7 0.0003043 ***
## treatment:hour               2    0.3014   0.3549      8     16 0.9295212    
## gender:hour                  1    0.2927   0.7243      4      7 0.6023742    
## treatment:gender:hour        2    0.5702   0.7976      8     16 0.6131884    
## phase:hour                   1    0.5496   0.4576      8      3 0.8324517    
## treatment:phase:hour         2    0.6637   0.2483     16      8 0.9914415    
## gender:phase:hour            1    0.6950   0.8547      8      3 0.6202076    
## treatment:gender:phase:hour  2    0.7928   0.3283     16      8 0.9723693    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

summary(av.ok, multivariate=FALSE)

## Univariate Type II Repeated-Measures ANOVA Assuming Sphericity
## 
##                                  SS num Df Error SS den Df       F    Pr(>F)
## treatment                   211.286      2  228.056     10  4.6323  0.037687
## gender                       58.286      1  228.056     10  2.5558  0.140974
## treatment:gender            130.241      2  228.056     10  2.8555  0.104469
## phase                       167.500      2   80.278     20 20.8651 1.274e-05
## treatment:phase              78.668      4   80.278     20  4.8997  0.006426
## gender:phase                  1.668      2   80.278     20  0.2078  0.814130
## treatment:gender:phase       10.221      4   80.278     20  0.6366  0.642369
## hour                        106.292      4   62.500     40 17.0067 3.191e-08
## treatment:hour                1.161      8   62.500     40  0.0929  0.999257
## gender:hour                   2.559      4   62.500     40  0.4094  0.800772
## treatment:gender:hour         7.755      8   62.500     40  0.6204  0.755484
## phase:hour                   11.083      8   96.167     80  1.1525  0.338317
## treatment:phase:hour          6.262     16   96.167     80  0.3256  0.992814
## gender:phase:hour             6.636      8   96.167     80  0.6900  0.699124
## treatment:gender:phase:hour  14.155     16   96.167     80  0.7359  0.749562
## 
## treatment                   *
## gender
## treatment:gender
## phase                       ***
## treatment:phase             **
## gender:phase
## treatment:gender:phase
## hour                        ***
## treatment:hour
## gender:hour
## treatment:gender:hour
## phase:hour
## treatment:phase:hour
## gender:phase:hour
## treatment:gender:phase:hour
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## 
## Mauchly Tests for Sphericity
## 
##                             Test statistic p-value
## phase                              0.74927 0.27282
## treatment:phase                    0.74927 0.27282
## gender:phase                       0.74927 0.27282
## treatment:gender:phase             0.74927 0.27282
## hour                               0.06607 0.00760
## treatment:hour                     0.06607 0.00760
## gender:hour                        0.06607 0.00760
## treatment:gender:hour              0.06607 0.00760
## phase:hour                         0.00478 0.44939
## treatment:phase:hour               0.00478 0.44939
## gender:phase:hour                  0.00478 0.44939
## treatment:gender:phase:hour        0.00478 0.44939
## 
## 
## Greenhouse-Geisser and Huynh-Feldt Corrections
##  for Departure from Sphericity
## 
##                              GG eps Pr(>F[GG])
## phase                       0.79953  7.323e-05 ***
## treatment:phase             0.79953    0.01223 *
## gender:phase                0.79953    0.76616
## treatment:gender:phase      0.79953    0.61162
## hour                        0.46028  8.741e-05 ***
## treatment:hour              0.46028    0.97879
## gender:hour                 0.46028    0.65346
## treatment:gender:hour       0.46028    0.64136
## phase:hour                  0.44950    0.34573
## treatment:phase:hour        0.44950    0.94019
## gender:phase:hour           0.44950    0.58903
## treatment:gender:phase:hour 0.44950    0.64634
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##                              HF eps Pr(>F[HF])
## phase                       0.92786  2.388e-05 ***
## treatment:phase             0.92786    0.00809 **
## gender:phase                0.92786    0.79845
## treatment:gender:phase      0.92786    0.63200
## hour                        0.55928  2.014e-05 ***
## treatment:hour              0.55928    0.98877
## gender:hour                 0.55928    0.69115
## treatment:gender:hour       0.55928    0.66930
## phase:hour                  0.73306    0.34405
## treatment:phase:hour        0.73306    0.98047
## gender:phase:hour           0.73306    0.65524
## treatment:gender:phase:hour 0.73306    0.70801
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

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

## Type II Repeated Measures MANOVA Tests: Roy test statistic
##               Df test stat approx F num Df den Df    Pr(>F)    
## measure        1    86.203  1293.04      2     30 < 2.2e-16 ***
## group:measure  2     0.356     5.52      2     31  0.008906 ** 
## month          1     9.407    65.85      4     28 7.807e-14 ***
## group:month    2     1.772    12.84      4     29 3.909e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

## mixed-effects models

## mixed-effects models examples:

\dontrun{
	library(nlme)
	example(lme)
	Anova(fm2)
}

## Analysis of Deviance Table (Type II tests)
## 
## Response: distance
##            Df    Chisq Pr(>Chisq)    
## age         1 114.8383  < 2.2e-16 ***
## Sex         1   9.2921   0.002301 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

\dontrun{
	library(lme4)
	example(lmer)
	Anova(gm1)
}

## Analysis of Deviance Table (Type II tests)
## 
## Response: cbind(incidence, size - incidence)
##           Df  Chisq Pr(>Chisq)    
## period     3 25.326  1.319e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

}
\keyword{htest}
\keyword{models}
\keyword{regression}
