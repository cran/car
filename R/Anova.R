#-------------------------------------------------------------------------------
# Revision history:
# 2009-01-05: bug fix in Anova_II_lm(). J. Fox
# 2009-01-16: Cox models with clusters now handled. J. Fox
# 2009-09-16: reworked glm and lm methods to handle aliased parameters. J. Fox
# 2009-09-30: renamed "Anova" to "Analysis of Deviance" in output for some methods. J. Fox
# 2009-12-22: modified Anova.mlm() to handle a user-supplied within-subject model matrix. J. Fox
# 2009-12-28: named the components of P in Anova_III_mlm(). John
# 2010-01-01: Anova_II_mlm() now hands off (again) to Anova_III_mlm() when there
#             is only an intercept in the between-subjects model
# 2010-02-17: Fixed bug that caused some models with aliased coefficients to fail. J. Fox
# 2010-06-14: added wcrossprod and allow use of observation weights in Anova.mlm()
# 2010-06-28: Fixed Anova() tables for coxph and survreg models 
#             (failed because of changes in survival package.
# 2011-01-21: Added functions for mixed models. J. Fox
# 2011-01-25: Fixed Anova.polr() and Anova.multinom() to work with models with only one term. J. Fox
# 2011-05-19: local fixef() to avoid nlme/lme4 issues. J. Fox
# 2011-05-11: changed order of columns in ANOVA tables for mixed models. J. Fox
# 2011-11-27: added Anova.svyglm(). J. Fox
# 2011-12-31: fixed bug in Anova.II(and III).F.glm() when na.exclude used. J. Fox
# 2012-02-28: added test.statistic argument to Anova.mer(). J.Fox
# 2012-03-02: fixed test abbreviation of test.statistic argument to Anova.default()
#             called by other Anova() methods. J. Fox
# 2013-06-17: modified summary.Anova.mlm(), introduced print.summary.Anova.mlm(),
#             adapting code contributed by Gabriel Baud-Bovy. J. Fox
# 2013-06-20: added Anova.merMod() method. J. Fox
# 2013-06-22: tweaks to local fixef(). J. Fox
# 2013-06-22: test argument uniformly uses "Chisq" rather than "chisq". J. Fox
# 2013-08-19: replaced calls to print.anova(). J. Fox
# 2014-08-17: added calls to requireNamespace() and :: where needed (doesn't work for pbkrtest). J. Fox
# 2014-08-18: fixed bugs in Anova.survreg() for types II, III LR tests and Wald tests. J. Fox
# 2014-09-23: added Anova.rlm(). J. Fox
# 2014-10-10: removed MASS:: from calls to polr(). John
# 2014-12-18: check that residual df and SS are nonzero in Anova.lm(). John
# 2015-01-27: vcovAdj() and methods now imported from pbkrtest. John
# 2015-02-18: force evaluation of vcov. when it's a function. John
# 2015-04-30: don't allow error.estimate="dispersion" for F-tests in binomial
#             and Poission GLMs. John
# 2015-08-29: fixed Anova() for coxph models with clusters. John
# 2015-09-04: added support for coxme models. John
# 2015-09-11: modified Anova.default() to work with vglm objects from VGAM. John
# 2015-09-15: fixed Anova.default() so that F-tests work again. John
# 2015-11-13: modify Anova.coxph() to take account of method/ties argument. John
# 2016-06-03: added SSP and SSPE args to print.summary.Anova.mlm(). John
# 2016-06-25: added code to optionally print univariate ANOVAs for a mlm. John
# 2017-02-16: replace polr() calls with MASS::polr(), multinom() with nnet::multinom(),
#             vcovAdj() with pbkrtest::vcovAdj(). John
# 2017-03-08: fixed bug in print.summary.Anova.mlm(). John
# 2017-11-07: added complete=FALSE to vcov() and vcov.() calls. John
# 2017-11-24: small improvements to output. John
# 2017-11-29: further fixed to vcov() and vcov.() calls. John
# 2018-01-15: Anova.multinom() now works with response matrix. JF
# 2018-02-11: If there are aliased coefs in lm object, treat as GLM. JF
# 2018-04-04: pass ... arguments through print() methods. Follows comments by Egor Katkov. JF
# 2019-10-16: modify Anova.coxph() and Anova.default()  for coxph() models with strata (or clusters)
#             (following problem reported by Susan Galloway Hilsenbeck). JF
# 2019-02-17: fix Anova.lme() to work with models without intercepts (to fix bug reported by Benjamin Tyner). JF
# 2020-04-01: fix Anova.coxph() to work with weights (to fix bug reported by Daniel Morillo Cuadrado)
# 2020-05-27: tweak to handling of Anova.coxph Wald tests. JF
# 2020-12-07: Standardize handling of vcov. arg
# 2020-12-18: fix Anova.lme() so that it handles missing factor levels. JF
# 2020-12-18: make assignVector() generic; add default and svyolr methods;
#             add unexported svyolr methods for coef() and vcov();
#             all this to make Anova() and linearHypothesis() work with svyolr. JF
# 2021-04-07: fix Anova.lm() so that SSs are computed when vcov. not specified. JF
# 2021-06-12: vcov. arg. now works for mer models.
# 2021-06-14: further fixes to vcov. arg for Anova.mer(). JF
#             introduced vcov. arg to Anova.glm(). JF
# 2021-06-16: Fix imatrix arg to Anova.mlm() (contribution of Benedikt Langenberg).JF
# 2021-06-19: make sure that calls to anova() for survival::survreg() models return "anova" objects. JF
# 2022-01-17,18: handle singularities better in Anova.mlm() (suggestion of Marius Barth)
# 2922-04-24: introduce new error.df argument for linearHypothesis.default(). JF
# 2022-06-07: Added Anova.svycoxph(). JF
# 2022-07-22: Fix bug in Anova.survreg() for Wald tests (reported by Megan Taylor Jones). JF
# 2022-07-22: Make Anova.lm() more robust when there are aliased coefficients (following report by Taiwo Fagbohungbe). JF
# 2022-07-27: Tweaked the last fix so the tolerance for deciding rank is the same for the lm model and the temporary glm model. SW
# 2023-10-03: Suppress LR tests for "coxph" models using the tt argument (following bug report by Ken Beath). JF
# 2024-05-08: Added Anova.clm() and Anova.clmm() methods (and supporting functions) (following report by Karl Ove Hufthammer). JF
# 2024-05-14: Rename internal functions to replace .s with _s. JF
# 2024-09-19: model.matrix.lme() now handles contrasts and xlev correctly, fixing a bug in Anova.lme() (reported by Ben Bolker). JF
#-------------------------------------------------------------------------------

# Type II and III tests for linear, generalized linear, and other models (J. Fox)

ConjComp <- function(X, Z = diag( nrow(X)), ip = diag(nrow(X))) {
  # This function by Georges Monette
  # finds the conjugate complement of the proj of X in span(Z) wrt
  #    inner product ip
  # - assumes Z is of full column rank
  # - projects X conjugately wrt ip into span Z
  xq <- qr(t(Z) %*% ip %*% X)
  if (xq$rank == 0) return(Z)
  Z %*% qr.Q(xq, complete = TRUE) [ ,-(1:xq$rank)] 
}

relatives <- function(term, names, factors){
  is.relative <- function(term1, term2) {
    all(!(factors[,term1]&(!factors[,term2])))
  }
  if(length(names) == 1) return(NULL)
  which.term <- which(term==names)
  (1:length(names))[-which.term][sapply(names[-which.term], 
                                        function(term2) is.relative(term, term2))]
}

lm2glm <- function(mod){
  Data <- getModelData(mod)
  wts <- weights(mod)
  Data$..wts.. <- if (is.null(wts)) rep(1, nrow(Data)) else wts
  form <- formula(mod)
  eps <- 1000 * (if(is.null(mod$call$tol)) 1e-7 else mod$call$tol)
  glm(form, weights=..wts.., data=Data, control=list(epsilon=eps))
}

globalVariables("..wts..")

Anova <- function(mod, ...){
  UseMethod("Anova", mod)
}

# linear models

Anova.lm <- function(mod, error, type=c("II","III", 2, 3), 
                     white.adjust=c(FALSE, TRUE, "hc3", "hc0", "hc1", "hc2", "hc4"),
                     vcov.=NULL, singular.ok, ...){
  if (!is.null(vcov.)) message("Coefficient covariances computed by ", deparse(substitute(vcov.)))
  if (!missing(white.adjust)) message("Coefficient covariances computed by hccm()")
  if (df.residual(mod) == 0) stop("residual df = 0")
  if (deviance(mod) < sqrt(.Machine$double.eps)) stop("residual sum of squares is 0 (within rounding error)")
  type <- as.character(type)
  white.adjust <- as.character(white.adjust)
  type <- match.arg(type)
  white.adjust <- match.arg(white.adjust)
  if (missing(singular.ok)){
    singular.ok <- type == "2" || type == "II"
  }
  if (has_intercept(mod) && length(coef(mod)) == 1 
      && (type == "2" || type == "II")) {
    type <- "III"
    warning("the model contains only an intercept: Type III test substituted")
  }
  if (any(is.na(coef(mod))) && singular.ok){
    if ((white.adjust != "FALSE") || (!is.null(vcov.)))
      stop("non-standard coefficient covariance matrix\n  may not be used for model with aliased coefficients")
    message("Note: model has aliased coefficients\n      sums of squares computed by model comparison")
    result <- Anova(lm2glm(mod), type=type, singular.ok=TRUE, test.statistic="F", ...)
    heading <- attributes(result)$heading
    if (type == "2") type <- "II"
    if (type == "3") type <- "III"
    attr(result, "heading") <- c(paste("Anova Table (Type", type, "tests)"), "", heading[2])
    return(result)
  }
  if (white.adjust != "FALSE"){
    if (white.adjust == "TRUE") white.adjust <- "hc3" 
    return(Anova.default(mod, type=type, vcov.=hccm(mod, type=white.adjust), test.statistic="F", 
                         singular.ok=singular.ok, ...))
  }
  else if (!is.null(vcov.)) return(Anova.default(mod, type=type, vcov.=vcov., test.statistic="F", 
                                                 singular.ok=singular.ok, ...))
  switch(type,
         II=Anova_II_lm(mod, error, singular.ok=singular.ok, ...),
         III=Anova_III_lm(mod, error, singular.ok=singular.ok, ...),
         "2"=Anova_II_lm(mod, error, singular.ok=singular.ok, ...),
         "3"=Anova_III_lm(mod, error, singular.ok=singular.ok,...))
}

Anova.aov <- function(mod, ...){
  class(mod) <- "lm"
  Anova.lm(mod, ...)
}

Anova_II_lm <- function(mod, error, singular.ok=TRUE, ...){
  if (!missing(error)){
    sumry <- summary(error, corr=FALSE)
    s2 <- sumry$sigma^2
    error.df <- error$df.residual
    error.SS <- s2*error.df
  }
  SS.term <- function(term){
    which.term <- which(term == names)
    subs.term <- which(assign == which.term)
    relatives <- relatives(term, names, fac)
    subs.relatives <- NULL
    for (relative in relatives) 
      subs.relatives <- c(subs.relatives, which(assign == relative))
    hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
    hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
    hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
    hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]
    hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
    else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov(mod, complete=FALSE)))
    hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
                                              function(x) all(x == 0)), , drop=FALSE]
    if (nrow(hyp.matrix.term) == 0)
      return(c(SS=NA, df=0))
    lh <- linearHypothesis(mod, hyp.matrix.term, 
                           singular.ok=singular.ok, ...)
    abs(c(SS=lh$"Sum of Sq"[2], df=lh$Df[2]))
  }
  not.aliased <- !is.na(coef(mod))
  if (!singular.ok && !all(not.aliased))
    stop("there are aliased coefficients in the model")
  fac <- attr(terms(mod), "factors")
  intercept <- has_intercept(mod)
  I.p <- diag(length(coefficients(mod)))
  assign <- mod$assign
  assign[!not.aliased] <- NA
  names <- term.names(mod)
  if (intercept) names <-names[-1]
  n.terms <- length(names)
  p <- df <- f <- SS <- rep(0, n.terms + 1)
  sumry <- summary(mod, corr = FALSE)
  SS[n.terms + 1] <- if (missing(error)) sumry$sigma^2*mod$df.residual 
  else error.SS   
  df[n.terms + 1] <- if (missing(error)) mod$df.residual else error.df
  p[n.terms + 1] <- f[n.terms + 1] <- NA
  for (i in 1:n.terms){
    ss <- SS.term(names[i])
    SS[i] <- ss["SS"]
    df[i] <- ss["df"]
    f[i] <- df[n.terms+1]*SS[i]/(df[i]*SS[n.terms + 1])
    p[i] <- pf(f[i], df[i], df[n.terms + 1], lower.tail = FALSE)
  }    
  result <- data.frame(SS, df, f, p)
  row.names(result) <- c(names,"Residuals")
  names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Anova Table (Type II tests)\n", 
                               paste("Response:", responseName(mod)))
  result
}

# type III

Anova_III_lm <- function(mod, error, singular.ok=FALSE, ...){
  if (!missing(error)){
    error.df <- df.residual(error)
    error.SS <- deviance(error)
  }
  else {
    error.df <- df.residual(mod)
    error.SS <- deviance(mod)
  }
  intercept <- has_intercept(mod)
  I.p <- diag(length(coefficients(mod)))
  Source <- term.names(mod)
  n.terms <- length(Source)
  p <- df <- f <- SS <- rep(0, n.terms + 1)
  assign <- mod$assign
  not.aliased <- !is.na(coef(mod))
  if (!singular.ok && !all(not.aliased))
    stop("there are aliased coefficients in the model")
  indices <- 1:n.terms
  for (term in indices){
    subs <- which(assign == term - intercept)
    hyp.matrix <- I.p[subs,,drop=FALSE]
    hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
    hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]
    if (nrow(hyp.matrix) == 0){
      SS[term] <- NA
      df[term] <- 0
      f[term] <- NA
      p[term] <- NA
    }
    else {
      test <- linearHypothesis(mod, hyp.matrix, singular.ok=singular.ok, ...)
      SS[term] <- test$"Sum of Sq"[2]
      df[term] <- test$"Df"[2]
    }
  }
  index.error <- n.terms + 1
  Source[index.error] <- "Residuals"
  SS[index.error] <- error.SS
  df[index.error] <- error.df
  f[indices] <- (SS[indices]/df[indices])/(error.SS/error.df)
  p[indices] <- pf(f[indices], df[indices], error.df, lower.tail=FALSE)
  p[index.error] <- f[index.error] <- NA
  result <- data.frame(SS, df, f, p)
  row.names(result) <- Source
  names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Anova Table (Type III tests)\n", paste("Response:", responseName(mod)))
  result
}

# generalized linear models

Anova.glm <- function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald", "F"),
                      error, error.estimate=c("pearson", "dispersion", "deviance"), 
                      vcov.=vcov(mod, complete=TRUE), singular.ok, ...){
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)
  error.estimate <- match.arg(error.estimate)
  if (!missing(vcov.)) {
    if (test.statistic != "Wald"){
      warning(paste0('test.statistic="', test.statistic,
                     '"; vcov. argument ignored'))
    } else {
      message("Coefficient covariances computed by ", deparse(substitute(vcov.)))
    }
  }
  vcov. <- getVcov(vcov., mod)
  if (has_intercept(mod) && length(coef(mod)) == 1 
      && (type == "2" || type == "II")) {
    type <- "III"
    warning("the model contains only an intercept: Type III test substituted")
  }
  if (missing(singular.ok)){
    singular.ok <- type == "2" || type == "II"
  }
  switch(type,
         II=switch(test.statistic,
                   LR=Anova_II_LR_glm(mod, singular.ok=singular.ok),
                   Wald=Anova.default(mod, type="II", singular.ok=singular.ok, vcov.=vcov.),
                   F=Anova_II_F_glm(mod, error, error.estimate, singular.ok=singular.ok)),
         III=switch(test.statistic,
                    LR=Anova_III_LR_glm(mod, singular.ok=singular.ok),
                    Wald=Anova.default(mod, type="III", singular.ok=singular.ok, vcov.=vcov.),
                    F=Anova_III_F_glm(mod, error, error.estimate, singular.ok=singular.ok)),
         "2"=switch(test.statistic,
                    LR=Anova_II_LR_glm(mod, singular.ok=singular.ok),
                    Wald=Anova.default(mod, type="II", singular.ok=singular.ok, vcov.=vcov.),
                    F=Anova_II_F_glm(mod, error, error.estimate, singular.ok=singular.ok)),
         "3"=switch(test.statistic,
                    LR=Anova_III_LR_glm(mod, singular.ok=singular.ok),
                    Wald=Anova.default(mod, type="III", singular.ok=singular.ok, vcov.=vcov.),
                    F=Anova_III_F_glm(mod, error, error.estimate, singular.ok=singular.ok)))
}


# type III

# LR test

Anova_III_LR_glm <- function(mod, singular.ok=FALSE, ...){
  if (!singular.ok && any(is.na(coef(mod))))
    stop("there are aliased coefficients in the model")
  Source <- if (has_intercept(mod)) term.names(mod)[-1]
  else term.names(mod)
  n.terms <- length(Source)
  p <- df <- LR <- rep(0, n.terms)
  dispersion <- summary(mod, corr = FALSE)$dispersion
  deviance <- deviance(mod)/dispersion
  for (term in 1:n.terms){
    mod.1 <- drop1(mod, scope=eval(parse(text=paste("~",Source[term]))))
    df[term] <- mod.1$Df[2]
    LR[term] <- if (df[term] == 0) NA else (mod.1$Deviance[2]/dispersion)-deviance
    p[term] <- pchisq(LR[term], df[term], lower.tail = FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- Source
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova","data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", paste("Response:", responseName(mod)))
  result
}

# F test

Anova_III_F_glm <- function(mod, error, error.estimate, singular.ok=FALSE, ...){
  if (!singular.ok && any(is.na(coef(mod))))
    stop("there are aliased coefficients in the model")
  fam <- family(mod)$family
  if ((fam == "binomial" || fam == "poisson") && error.estimate == "dispersion"){
    warning("dispersion parameter estimated from the Pearson residuals, not taken as 1")
    error.estimate <- "pearson"
  }
  if (missing(error)) error <- mod
  df.res <- df.residual(error)
  error.SS <- switch(error.estimate,
                     pearson=sum(residuals(error, "pearson")^2, na.rm=TRUE),
                     dispersion=df.res*summary(error, corr = FALSE)$dispersion,
                     deviance=deviance(error))
  Source <- if (has_intercept(mod)) term.names(mod)[-1]
  else term.names(mod)
  n.terms <- length(Source)
  p <- df <- f <- SS <-rep(0, n.terms+1)
  f[n.terms+1] <- p[n.terms+1] <- NA
  df[n.terms+1] <- df.res
  SS[n.terms+1] <- error.SS
  dispersion <- error.SS/df.res
  deviance <- deviance(mod)
  for (term in 1:n.terms){
    mod.1 <- drop1(mod, scope=eval(parse(text=paste("~",Source[term]))))
    df[term] <- mod.1$Df[2]
    SS[term] <- mod.1$Deviance[2] - deviance
    f[term] <- if (df[term] == 0) NA else (SS[term]/df[term])/dispersion
    p[term] <- pf(f[term], df[term], df.res, lower.tail = FALSE)
  }
  result <- data.frame(SS, df, f, p)
  row.names(result) <- c(Source, "Residuals")
  names(result) <- c("Sum Sq", "Df", "F values", "Pr(>F)")
  class(result) <- c("anova","data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", 
                               paste("Response:", responseName(mod)), 
                               paste("Error estimate based on",
                                     switch(error.estimate, 
                                            pearson="Pearson residuals", dispersion="estimated dispersion", 
                                            deviance="deviance"), "\n"))
  result
}

# type II

# LR test

Anova_II_LR_glm <- function(mod, singular.ok=TRUE, ...){
  if (!singular.ok && any(is.na(coef(mod))))
    stop("there are aliased coefficients in the model")
  # (some code adapted from drop1.glm)
  which.nms <- function(name) which(asgn == which(names == name))
  fac <- attr(terms(mod), "factors")
  names <- if (has_intercept(mod)) term.names(mod)[-1]
  else term.names(mod)
  n.terms <- length(names)
  X <- model.matrix(mod)
  y <- mod$y
  if (is.null(y)) y <- model.response(model.frame(mod), "numeric")
  wt <- mod$prior.weights
  if (is.null(wt)) wt <- rep(1, length(y))
  asgn <- attr(X, 'assign')
  df <- p <- LR <- rep(0, n.terms)
  dispersion <- summary(mod, corr = FALSE)$dispersion
  for (term in 1:n.terms){
    rels <- names[relatives(names[term], names, fac)]
    exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
    mod.1 <- glm.fit(X[, -exclude.1, drop = FALSE], y, wt, offset = mod$offset, 
                     family = mod$family, control = mod$control)
    dev.1 <- deviance(mod.1)
    mod.2 <- if (length(rels) == 0) mod
    else {
      exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
      glm.fit(X[, -exclude.2, drop = FALSE], y, wt, offset = mod$offset, 
              family = mod$family, control = mod$control)
    }
    dev.2 <- deviance(mod.2)
    df[term] <- df.residual(mod.1) - df.residual(mod.2)
    if (df[term] == 0) LR[term] <- p[term] <- NA
    else {
      LR[term] <- (dev.1 - dev.2)/dispersion
      p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
    }
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- 
    c("Analysis of Deviance Table (Type II tests)\n", paste("Response:", responseName(mod)))
  result
}


# F test

Anova_II_F_glm <- function(mod, error, error.estimate, singular.ok=TRUE, ...){
  # (some code adapted from drop1.glm)
  if (!singular.ok && any(is.na(coef(mod))))
    stop("there are aliased coefficients in the model")
  fam <- family(mod)$family
  if ((fam == "binomial" || fam == "poisson") && error.estimate == "dispersion"){
    warning("dispersion parameter estimated from the Pearson residuals, not taken as 1")
    error.estimate <- "pearson"
  }
  which.nms <- function(name) which(asgn == which(names == name))
  if (missing(error)) error <- mod
  df.res <- df.residual(error)
  error.SS <- switch(error.estimate,
                     pearson = sum(residuals(error, "pearson")^2, na.rm=TRUE),
                     dispersion = df.res*summary(error, corr = FALSE)$dispersion,
                     deviance = deviance(error))
  fac <- attr(terms(mod), "factors")
  names <- if (has_intercept(mod)) term.names(mod)[-1]
  else term.names(mod)
  n.terms <- length(names)
  X <- model.matrix(mod)
  y <- mod$y
  if (is.null(y)) y <- model.response(model.frame(mod), "numeric")
  wt <- mod$prior.weights
  if (is.null(wt)) wt <- rep(1, length(y))
  asgn <- attr(X, 'assign')
  p <- df <- f <- SS <- rep(0, n.terms+1)
  f[n.terms+1] <- p[n.terms+1] <- NA
  df[n.terms+1] <- df.res
  SS[n.terms+1] <- error.SS
  dispersion <- error.SS/df.res
  for (term in 1:n.terms){
    rels <- names[relatives(names[term], names, fac)]
    exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
    mod.1 <- glm.fit(X[, -exclude.1, drop = FALSE], y, wt, offset = mod$offset, 
                     family = mod$family, control = mod$control)
    dev.1 <- deviance(mod.1)
    mod.2 <- if (length(rels) == 0) mod
    else {
      exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
      glm.fit(X[, -exclude.2, drop = FALSE], y, wt, offset = mod$offset, 
              family = mod$family, control = mod$control)
    }
    dev.2 <- deviance(mod.2)
    df[term] <- df.residual(mod.1) - df.residual(mod.2)
    if (df[term] == 0) SS[term] <- f[term] <- p[term] <- NA
    else {
      SS[term] <- dev.1 - dev.2
      f[term] <- SS[term]/(dispersion*df[term])
      p[term] <- pf(f[term], df[term], df.res, lower.tail=FALSE)
    }
  }
  result <- data.frame(SS, df, f, p)
  row.names(result) <- c(names, "Residuals")
  names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
                               paste("Response:", responseName(mod)),
                               paste("Error estimate based on",
                                     switch(error.estimate,
                                            pearson="Pearson residuals", 
                                            dispersion="estimated dispersion", 
                                            deviance="deviance"), "\n"))
  result
}

# multinomial logit models (via multinom in the nnet package)

Anova.multinom <-
  function (mod, type = c("II", "III", 2, 3), ...)
  {
    type <- as.character(type)
    type <- match.arg(type)
    if (has_intercept(mod) && length(coef(mod)) == 1 
        && (type == "2" || type == "II")) {
      type <- "III"
      warning("the model contains only an intercept: Type III test substituted")
    }
    switch(type,
           II = Anova_II_multinom(mod, ...),
           III = Anova_III_multinom(mod, ...),
           "2" = Anova_II_multinom(mod, ...),
           "3" = Anova_III_multinom(mod, ...))
  }

Anova_II_multinom <- function (mod, ...)
{
  which.nms <- function(name) which(asgn == which(names ==
                                                    name))
  fac <- attr(terms(mod), "factors")
  names <- if (has_intercept(mod)) term.names(mod)[-1]
  else term.names(mod)
  n.terms <- length(names)
  X <- model.matrix(mod)
  y <- model.response(model.frame(mod))
  wt <- if (is.matrix(y)) rep(1, nrow(y)) else mod$weights
  asgn <- attr(X, "assign")
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  for (term in 1:n.terms) {
    rels <- names[relatives(names[term], names, fac)]
    exclude.1 <- as.vector(unlist(sapply(c(names[term], rels),
                                         which.nms)))
    mod.1 <-if (n.terms > 1) nnet::multinom(y ~ X[, -c(1, exclude.1)], weights=wt, trace=FALSE)
    else nnet::multinom(y ~ 1, weights=wt, race=FALSE)
    dev.1 <- deviance(mod.1)
    mod.2 <- if (length(rels) == 0)
      mod
    else {
      exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
      nnet::multinom(y ~ X[, -c(1, exclude.2)], weights=wt, trace=FALSE)
    }
    dev.2 <- deviance(mod.2)
    LR[term] <- dev.1 - dev.2
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n",
                               paste("Response:", responseName(mod)))
  result
}

Anova_III_multinom <- function (mod, ...)
{
  names <- if (has_intercept(mod)) term.names(mod)[-1]
  else term.names(mod)
  n.terms <- length(names)
  X <- model.matrix(mod)
  y <- model.response(model.frame(mod))
  wt <- if (is.matrix(y)) rep(1, nrow(y)) else mod$weights
  asgn <- attr(X, "assign")
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  deviance <- deviance(mod)
  for (term in 1:n.terms) {
    mod.1 <- if (n.terms > 1) nnet::multinom(y ~ X[, term != asgn][, -1], weights=wt, trace=FALSE)
    else nnet::multinom(y ~ 1, weights=wt, trace=FALSE)
    LR[term] <- deviance(mod.1) - deviance
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n",
                               paste("Response:", responseName(mod)))
  result
}


# proportional-odds logit models (via polr in the MASS package)

Anova.polr <- function (mod, type = c("II", "III", 2, 3), ...)
{
  type <- as.character(type)
  type <- match.arg(type)
  if (has_intercept(mod) && length(coef(mod)) == 1 
      && (type == "2" || type == "II")) {
    type <- "III"
    warning("the model contains only an intercept: Type III test substituted")
  }
  switch(type,
         II = Anova_II_polr(mod, ...),
         III = Anova_III_polr(mod, ...),
         "2" = Anova_II_polr(mod, ...),
         "3" = Anova_III_polr(mod, ...))
}

Anova_II_polr <- function (mod, ...)
{
  if (!requireNamespace("MASS")) stop("MASS package is missing")
  which.nms <- function(name) which(asgn == which(names ==
                                                    name))
  fac <- attr(terms(mod), "factors")
  names <- term.names(mod)
  n.terms <- length(names)
  X <- model.matrix(mod)
  y <- model.response(model.frame(mod))
  wt <- model.weights(model.frame(mod))
  asgn <- attr(X, "assign")
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  for (term in 1:n.terms) {
    rels <- names[relatives(names[term], names, fac)]
    exclude.1 <- as.vector(unlist(sapply(c(names[term], rels),
                                         which.nms)))
    mod.1 <- if (n.terms > 1) MASS::polr(y ~ X[, -c(1, exclude.1)], weights=wt)
    else MASS::polr(y ~ 1, weights=wt)
    dev.1 <- deviance(mod.1)
    mod.2 <- if (length(rels) == 0)
      mod
    else {
      exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
      MASS::polr(y ~ X[, -c(1, exclude.2)], weights=wt)
    }
    dev.2 <- deviance(mod.2)
    LR[term] <- dev.1 - dev.2
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n",
                               paste("Response:", responseName(mod)))
  result
}

Anova_III_polr <- function (mod, ...)
{
  if (!requireNamespace("MASS")) stop("MASS package is missing")
  names <- term.names(mod)
  n.terms <- length(names)
  X <- model.matrix(mod)
  y <- model.response(model.frame(mod))
  wt <- model.weights(model.frame(mod))
  asgn <- attr(X, "assign")
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  deviance <- deviance(mod)
  for (term in 1:n.terms) {
    mod.1 <- if (n.terms > 1) MASS::polr(y ~ X[, term != asgn][, -1], weights=wt)
    else MASS::polr(y ~ 1, weights=wt)
    LR[term] <- deviance(mod.1) - deviance
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n",
                               paste("Response:", responseName(mod)))
  result
}

# multivariate linear models

# the following 3 functions copied from the stats package (not exported from stats)

Pillai <- function (eig, q, df.res) {
  test <- sum(eig/(1 + eig))
  p <- length(eig)
  s <- min(p, q)
  n <- 0.5 * (df.res - p - 1)
  m <- 0.5 * (abs(p - q) - 1)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * n + s + 1
  c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s * tmp2)
}

Wilks <- function (eig, q, df.res) {
  test <- prod(1/(1 + eig))
  p <- length(eig)
  tmp1 <- df.res - 0.5 * (p - q + 1)
  tmp2 <- (p * q - 2)/4
  tmp3 <- p^2 + q^2 - 5
  tmp3 <- if (tmp3 > 0) 
    sqrt(((p * q)^2 - 4)/tmp3)
  else 1
  c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q, 
    p * q, tmp1 * tmp3 - 2 * tmp2)
}

HL <- function (eig, q, df.res) {
  test <- sum(eig)
  p <- length(eig)
  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (df.res - p - 1)
  s <- min(p, q)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * (s * n + 1)
  c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
}

Roy <- function (eig, q, df.res) {
  p <- length(eig)
  test <- max(eig)
  tmp1 <- max(p, q)
  tmp2 <- df.res - tmp1 + q
  c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
}

has_intercept.mlm <- function (model, ...) 
  any(row.names(coefficients(model)) == "(Intercept)")


Anova.mlm <- function(mod, type=c("II","III", 2, 3), SSPE, error.df, idata, 
                      idesign, icontrasts=c("contr.sum", "contr.poly"), imatrix,
                      test.statistic=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),...){
  if (any(is.na(coef(mod)))) 
    stop(if(!missing(idata)) "between-subjects ", "model is singular")
  wts <- if (!is.null(mod$weights)) mod$weights else rep(1, nrow(model.matrix(mod)))
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)
  if (missing(SSPE)) SSPE <- wcrossprod(residuals(mod), w=wts)
  if (missing(idata)) {
    idata <- NULL
    idesign <- NULL
  }
  if (missing(imatrix)) imatrix <- NULL
  error.df <- if (missing(error.df)) df.residual(mod)
  else error.df
  switch(type,
         II=Anova_II_mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
         III=Anova_III_mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
         "2"=Anova_II_mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...),
         "3"=Anova_III_mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test.statistic, ...))
}

Anova_III_mlm <- function(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...){
  intercept <- has_intercept(mod)
  p <- nrow(coefficients(mod))
  I.p <- diag(p)
  terms <- term.names(mod)
  n.terms <- length(terms)
  assign <- mod$assign
  if (is.null(idata) && is.null(imatrix)){
    if ((n.terms == 0) && intercept) {
      Test <- linearHypothesis(mod, 1, SSPE=SSPE, ...)
      result <- list(SSP=Test$SSPH, SSPE=SSPE, df=1, error.df=error.df,
                     terms="(Intercept)", repeated=FALSE, type="III", test=test)
      class(result) <- "Anova.mlm"
      return(result)
    }
    SSP <- as.list(rep(0, n.terms))
    df <- rep(0, n.terms)
    names(df) <- names(SSP) <- terms
    for (term in 1:n.terms){
      subs <- which(assign == term - intercept)
      hyp.matrix <- I.p[subs,,drop=FALSE]
      Test <- linearHypothesis(mod, hyp.matrix, SSPE=SSPE, ...)
      SSP[[term]] <- Test$SSPH
      df[term]<- length(subs)
    }
    result <- list(SSP=SSP, SSPE=SSPE, df=df, error.df=error.df, terms=terms,
                   repeated=FALSE, type="III", test=test)
  }
  else {
    if (!is.null(imatrix)){
      X.design <- do.call(cbind, imatrix)
      ncols <- sapply(imatrix, ncol)
      end <- cumsum(ncols)
      start <- c(1, (end + 1))[-(length(end) + 1)]
      cols <- mapply(seq, from=start, to=end)
      iterms <- names(end)
      names(cols) <- iterms
      itrms <- unlist(sapply(1:length(imatrix), function(x) replicate(ncol(imatrix[[x]]), x-1)))
      check.imatrix(X.design, itrms)
    }
    else {
      if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
      for (i in 1:length(idata)){
        if (is.null(attr(idata[,i], "contrasts"))){
          contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
          else icontrasts[1]
        }
      }
      X.design <- model.matrix(idesign, data=idata)
      i.intercept <- has_intercept(X.design)
      iterms <- term.names(idesign)
      if (i.intercept) iterms <- c("(Intercept)", iterms)
      check.imatrix(X.design)
    }
    n.tests <- n.terms*length(iterms)
    df <- rep(0, n.tests)
    singular <- rep(FALSE, n.tests)
    hnames <- rep("", n.tests)
    P <- SSPEH <- SSP <- as.list(df)
    i <- 0
    for (iterm in iterms){
      for (term in 1:n.terms){
        subs <- which(assign == term - intercept)
        hyp.matrix <- I.p[subs,,drop=FALSE]
        i <- i + 1
        Test <- linearHypothesis(mod, hyp.matrix, SSPE=SSPE, 
                                 idata=idata, idesign=idesign, icontrasts=icontrasts, iterms=iterm, 
                                 check.imatrix=FALSE, P=imatrix[[iterm]], singular.ok=TRUE, ...)
        SSP[[i]] <- Test$SSPH
        SSPEH[[i]] <- Test$SSPE
        P[[i]] <- Test$P
        df[i] <- length(subs)
        hnames[i] <- if (iterm == "(Intercept)") terms[term]
        else if (terms[term] == "(Intercept)") iterm
        else paste(terms[term], ":", iterm, sep="")
        singular[i] <- Test$singular
      }
    }
    names(singular) <- names(df) <- names(SSP) <- names(SSPEH) <- names(P) <- hnames
    result <- list(SSP=SSP, SSPE=SSPEH, P=P, df=df, error.df=error.df,
                   terms=hnames, repeated=TRUE, type="III", test=test, 
                   idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix,
                   singular=singular)       
  }
  class(result) <- "Anova.mlm"
  result
}

Anova_II_mlm <- function(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...){
  wts <- if (!is.null(mod$weights)) mod$weights else rep(1, nrow(model.matrix(mod)))
  V <- solve(wcrossprod(model.matrix(mod), w=wts))
  SSP.term <- function(term, iterm){
    which.term <- which(term == terms)
    subs.term <- which(assign == which.term)
    relatives <- relatives(term, terms, fac)
    subs.relatives <- NULL
    for (relative in relatives) subs.relatives <- c(subs.relatives, which(assign==relative))
    hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
    hyp.matrix.2 <- I.p[c(subs.relatives, subs.term),,drop=FALSE]
    if (missing(iterm)){
      SSP1 <- if (length(subs.relatives) == 0) 0 
      else linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, singular.ok=TRUE, ...)$SSPH
      SSP2 <- linearHypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, singular.ok=TRUE, ...)$SSPH
      return(SSP2 - SSP1)
    }
    else {
      SSP1 <- if (length(subs.relatives) == 0) 0 
      else linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
                            idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, P=imatrix[[iterm]], singular.ok=TRUE, ...)$SSPH
      lh2 <- linearHypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, 
                              idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, P=imatrix[[iterm]], singular.ok=TRUE, ...)
      return(list(SSP = lh2$SSPH - SSP1, SSPE=lh2$SSPE, P=lh2$P, singular=lh2$singular))
    }
  }
  fac <- attr(terms(mod), "factors")
  intercept <- has_intercept(mod)
  p <- nrow(coefficients(mod))
  I.p <- diag(p)
  assign <- mod$assign
  terms <- term.names(mod)
  if (intercept) terms <- terms[-1]
  n.terms <- length(terms)
  if (n.terms == 0){
    message("Note: model has only an intercept; equivalent type-III tests substituted.")
    return(Anova_III_mlm(mod, SSPE, error.df, idata, idesign, icontrasts, imatrix, test, ...))
  }
  if (is.null(idata) && is.null(imatrix)){
    SSP <- as.list(rep(0, n.terms))
    df <- rep(0, n.terms)
    names(df) <- names(SSP) <- terms
    for (i in 1:n.terms){
      SSP[[i]] <- SSP.term(terms[i])
      df[i]<- df.terms(mod, terms[i])
    }    
    result <- list(SSP=SSP, SSPE=SSPE, df=df, error.df=error.df, terms=terms,
                   repeated=FALSE, type="II", test=test)
  }
  else {
    if (!is.null(imatrix)){
      X.design <- do.call(cbind, imatrix)
      ncols <- sapply(imatrix, ncol)
      end <- cumsum(ncols)
      start <- c(1, (end + 1))[-(length(end) + 1)]
      cols <- mapply(seq, from=start, to=end)
      iterms <- names(end)
      names(cols) <- iterms
      itrms <- unlist(sapply(1:length(imatrix), function(x) replicate(ncol(imatrix[[x]]), x-1)))
      check.imatrix(X.design, itrms)     
    }
    else {
      if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
      for (i in 1:length(idata)){
        if (is.null(attr(idata[,i], "contrasts"))){
          contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
          else icontrasts[1]
        }
      }
      X.design <- model.matrix(idesign, data=idata)
      iintercept <- has_intercept(X.design)
      iterms <- term.names(idesign)
      if (iintercept) iterms <- c("(Intercept)", iterms)
      check.imatrix(X.design)
    }
    n.tests <-(n.terms + intercept)*length(iterms)
    df <- rep(0, n.tests)
    singular <- rep(FALSE, n.tests)
    hnames <- rep("", length(df))
    P <- SSPEH <- SSP <- as.list(df)
    i <- 0
    for (iterm in iterms){
      if (intercept){
        i <- i + 1
        hyp.matrix.1 <- I.p[-1,,drop=FALSE]
        SSP1 <- linearHypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
                                 idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, 
                                 check.imatrix=FALSE, P=imatrix[[iterm]], singular.ok=TRUE, ...)$SSPH
        lh2 <- linearHypothesis(mod, I.p, SSPE=SSPE, V=V, 
                                idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, 
                                check.imatrix=FALSE, P=imatrix[[iterm]], singular.ok=TRUE, ...)
        SSP[[i]] <- lh2$SSPH - SSP1
        SSPEH[[i]] <- lh2$SSPE
        P[[i]] <- lh2$P
        singular[i] <- lh2$singular
        df[i] <- 1
        hnames[i] <- iterm
      }
      for (term in 1:n.terms){
        subs <- which(assign == term)
        i <- i + 1
        Test <- SSP.term(terms[term], iterm)
        SSP[[i]] <- Test$SSP
        SSPEH[[i]] <- Test$SSPE
        P[[i]] <- Test$P
        singular[i] <- Test$singular
        df[i]<- length(subs)
        hnames[i] <- if (iterm == "(Intercept)") terms[term]
        else paste(terms[term], ":", iterm, sep="")
      }
    }
    names(singular) <- names(df) <- names(P) <- names(SSP) <- names(SSPEH) <- hnames
    result <- list(SSP=SSP, SSPE=SSPEH, P=P, df=df, error.df=error.df,
                   terms=hnames, repeated=TRUE, type="II", test=test,
                   idata=idata, idesign=idesign, icontrasts=icontrasts, imatrix=imatrix,
                   singular=singular)       
  }
  class(result) <- "Anova.mlm"
  result
}

print.Anova.mlm <- function(x, ...){
  if ((!is.null(x$singular)) && any(x$singular)) stop("singular error SSP matrix; multivariate tests unavailable\ntry summary(object, multivariate=FALSE)")
  test <- x$test
  repeated <- x$repeated
  ntests <- length(x$terms)
  tests <- matrix(NA, ntests, 4)
  if (!repeated) SSPE.qr <- qr(x$SSPE) 
  for (term in 1:ntests){
    # some of the code here adapted from stats:::summary.manova
    eigs <- Re(eigen(qr.coef(if (repeated) qr(x$SSPE[[term]]) else SSPE.qr,
                             x$SSP[[term]]), symmetric = FALSE)$values)
    tests[term, 1:4] <- switch(test,
                               Pillai = Pillai(eigs, x$df[term], x$error.df),
                               Wilks = Wilks(eigs, x$df[term], x$error.df),
                               "Hotelling-Lawley" = HL(eigs, x$df[term], x$error.df),
                               Roy = Roy(eigs, x$df[term], x$error.df))
  }
  ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
  ok <- !is.na(ok) & ok
  tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3], tests[ok, 4], 
                                 lower.tail = FALSE))
  rownames(tests) <- x$terms
  colnames(tests) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)")
  tests <- structure(as.data.frame(tests), 
                     heading = paste("\nType ", x$type, if (repeated) " Repeated Measures",
                                     " MANOVA Tests: ", test, " test statistic", sep=""), 
                     class = c("anova", "data.frame"))
  print(tests, ...)      
  invisible(x)
}


# summary.Anova.mlm and print.summary.Anova.mlm methods
#  with contributions from Gabriel Baud-Bovy
summary.Anova.mlm <- function (object, test.statistic, univariate=object$repeated, multivariate=TRUE, p.adjust.method, ...) {
  GG <- function(SSPE, P) { # Greenhouse-Geisser correction
    p <- nrow(SSPE)
    if (p < 2) 
      return(NA)
    lambda <- eigen(SSPE %*% solve(t(P) %*% P))$values
    lambda <- lambda[lambda > 0]
    ((sum(lambda)/p)^2)/(sum(lambda^2)/p)
  }
  HF <- function(gg, error.df, p) { # Huynh-Feldt correction
    ((error.df + 1) * p * gg - 2)/(p * (error.df - p * gg))
  }
  mauchly <- function(SSD, P, df) {
    # most of this function borrowed from stats:::mauchly.test.SSD
    if (nrow(SSD) < 2) 
      return(c(NA, NA))
    Tr <- function(X) sum(diag(X))
    p <- nrow(P)
    I <- diag(p)
    Psi <- t(P) %*% I %*% P
    B <- SSD
    pp <- nrow(SSD)
    U <- solve(Psi, B)
    n <- df
    logW <- log(det(U)) - pp * log(Tr(U/pp))
    rho <- 1 - (2 * pp^2 + pp + 2)/(6 * pp * n)
    w2 <- (pp + 2) * (pp - 1) * (pp - 2) * (2 * pp^3 + 6 * 
                                              pp^2 + 3 * p + 2)/(288 * (n * pp * rho)^2)
    z <- -n * rho * logW
    f <- pp * (pp + 1)/2 - 1
    Pr1 <- pchisq(z, f, lower.tail = FALSE)
    Pr2 <- pchisq(z, f + 4, lower.tail = FALSE)
    pval <- Pr1 + w2 * (Pr2 - Pr1)
    c(statistic = c(W = exp(logW)), p.value = pval)
  }
  if (missing(test.statistic)) 
    test.statistic <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
  test.statistic <- match.arg(test.statistic, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"), several.ok = TRUE)
  nterms <- length(object$terms)
  summary.object <- list(type=object$type, repeated=object$repeated, 
                         multivariate.tests=NULL, univariate.tests=NULL, 
                         pval.adjustments=NULL, sphericity.tests=NULL)
  if (multivariate){
    summary.object$multivariate.tests <- vector(nterms, mode="list")
    names(summary.object$multivariate.tests) <- object$terms 
    summary.object$SSPE <- object$SSPE
    for (term in 1:nterms) {
      hyp <- list(SSPH = object$SSP[[term]], 
                  SSPE = if (object$repeated) object$SSPE[[term]] else object$SSPE, 
                  P = if (object$repeated) object$P[[term]] else NULL, 
                  test = test.statistic, df = object$df[term], 
                  df.residual = object$error.df, title = object$terms[term],
                  singular = if (!is.null(object$singular)) object$singular[term]
                  )
      class(hyp) <- "linearHypothesis.mlm"
      summary.object$multivariate.tests[[term]] <- hyp
    }
  }
  if (object$repeated && univariate) {
    singular <- object$singular
    error.df <- object$error.df
    table <- matrix(0, nterms, 6)
    table2 <- matrix(0, nterms, 4)
    table3 <- matrix(0, nterms, 2)
    rownames(table3) <- rownames(table2) <- rownames(table) <- object$terms
    colnames(table) <- c("Sum Sq", "num Df", "Error SS", "den Df", "F value", "Pr(>F)")
    colnames(table2) <- c("GG eps", "Pr(>F[GG])", "HF eps","Pr(>F[HF])")
    colnames(table3) <- c("Test statistic", "p-value")
    if (any(singular))
      warning("one or more error SSP matrix:\ncorresponding non-sphericity tests and corrections not available")
    for (term in 1:nterms) {
      SSP <- object$SSP[[term]]
      SSPE <- object$SSPE[[term]]
      P <- object$P[[term]]
      p <- ncol(P)
      PtPinv <- solve(t(P) %*% P)
      gg <- if (!singular[term]) GG(SSPE, P) else NA
      table[term, "Sum Sq"] <- sum(diag(SSP %*% PtPinv))
      table[term, "Error SS"] <- sum(diag(SSPE %*% PtPinv))
      table[term, "num Df"] <- object$df[term] * p
      table[term, "den Df"] <- error.df * p
      table[term, "F value"] <- (table[term, "Sum Sq"]/table[term, "num Df"])/
        (table[term, "Error SS"]/table[term, "den Df"])
      table[term, "Pr(>F)"] <- pf(table[term, "F value"], table[term, "num Df"], table[term, "den Df"], 
                                  lower.tail = FALSE)
      table2[term, "GG eps"] <- gg
      table2[term, "HF eps"] <- if (!singular[term]) HF(gg, error.df, p) else NA
      table3[term, ] <- if (!singular[term]) mauchly(SSPE, P, object$error.df) else NA
    }
    table3 <- na.omit(table3)
    if (nrow(table3) > 0) {
      table2[, "Pr(>F[GG])"] <- pf(table[, "F value"], table2[, "GG eps"] * 
                                     table[, "num Df"], table2[, "GG eps"] * table[, "den Df"], 
                                   lower.tail = FALSE)
      table2[, "Pr(>F[HF])"] <- pf(table[, "F value"], pmin(1, table2[, "HF eps"]) * 
                                     table[, "num Df"], pmin(1, table2[, "HF eps"]) * table[, "den Df"], 
                                   lower.tail = FALSE)
      table2 <- na.omit(table2)
      if (any(table2[, "HF eps"] > 1)) warning("HF eps > 1 treated as 1")
    }
    class(table3) <- class(table) <- "anova"
    summary.object$univariate.tests <- table
    summary.object$pval.adjustments <- table2
    summary.object$sphericity.tests <- table3
  }
  if (!object$repeated && univariate) {
    SS <- sapply(object$SSP, diag)
    SSE <- diag(object$SSPE)
    df <- object$df
    dfe <- object$error.df
    F <- (SS/df)/(SSE/dfe)
    SS <- cbind(SS, residuals=SSE)
    SS <- rbind(df=c(df, residuals=dfe), SS)
    p <- pf(F, df, dfe, lower.tail=FALSE)
    result <- list(SS=t(SS), F=t(F), p=t(p), type=object$type)
    if (!missing(p.adjust.method)){
      if (isTRUE(p.adjust.method)) p.adjust.method <- "holm"
      p.adj <- apply(p, 2, p.adjust, method=p.adjust.method)
      result$p.adjust <- t(p.adj)
      result$p.adjust.method <- p.adjust.method
    }
    class(result) = "univaov"
    summary.object$univaov <- result
  }
  class(summary.object) <- "summary.Anova.mlm"
  summary.object
}

print.summary.Anova.mlm <- function(x, digits = getOption("digits"), SSP=TRUE, SSPE=SSP, ... ) {
  if (!is.null(x$multivariate.tests)) {
    cat(paste("\nType ", x$type, if (x$repeated) 
      " Repeated Measures", " MANOVA Tests:\n", sep = ""))
    if ((!x$repeated) && SSPE) {
      cat("\nSum of squares and products for error:\n")
      print(x$SSPE, digits = digits, ...)
    }
    for (term in 1:length(x$multivariate.tests)) {
      cat(paste("\n------------------------------------------\n", 
                "\nTerm:", names(x$multivariate.tests)[term], "\n"))
      print(x$multivariate.tests[[term]], digits = digits, SSP=SSP, SSPE=FALSE, ...)
    }
  }
  if  (!is.null(x$univariate.tests)) {
    cat("\nUnivariate Type", x$type, "Repeated-Measures ANOVA Assuming Sphericity\n\n")
    print(x$univariate.tests)
    if (nrow(x$sphericity.tests) > 0) {
      cat("\n\nMauchly Tests for Sphericity\n\n")
      print(x$sphericity.tests)
      cat("\n\nGreenhouse-Geisser and Huynh-Feldt Corrections\n", 
          "for Departure from Sphericity\n\n")
      table <- x$pval.adjustments[, 1:2, drop = FALSE]
      class(table) <- "anova"
      print(table, ...)
      cat("\n")
      table <- x$pval.adjustments[, 3:4, drop = FALSE]
      class(table)
      print(table, ...)
    }
  }
  if (!is.null(x$univaov)){
    print(x$univaov, ...)
  }
  invisible(x)
}

print.univaov <- function(x, digits = max(getOption("digits") - 2L, 3L), 
                          style=c("wide", "long"), 
                          by=c("response", "term"),
                          ...){
  style <- match.arg(style)
  if (style == "wide") {
    cat("\n Type", x$type, "Sums of Squares\n")
    print(x$SS, digits=digits)
    cat("\n F-tests\n")
    F <- x$F
    print(round(F, 2))
    cat("\n p-values\n")
    p <- format.pval(x$p)
    p <- matrix(p, nrow=nrow(F))
    rownames(p) <- rownames(F)
    colnames(p) <- colnames(F)
    print(p, quote=FALSE)
    if (!is.null(x$p.adjust)){
      cat("\n p-values adjusted (by term) for simultaneous inference by", x$p.adjust.method, "method\n")
      p.adjust <- format.pval(x$p.adjust)
      p.adjust <- matrix(p.adjust, nrow=nrow(F))
      rownames(p.adjust) <- rownames(F)
      colnames(p.adjust) <- colnames(F)
      print(p.adjust, quote=FALSE)
    }
  }
  else {
    x.df <- as.data.frame(x, by=by)
    x.df$F <- round(x.df$F, 2)
    x.df$p <- format.pval(x.df$p)
    if (!is.null(x$p.adjust)) x.df$"adjusted p" <- format.pval(x.df$"adjusted p")
    cat("\n Type", x$type, "Sums of Squares and F tests\n")
    print(x.df, quote=FALSE, digits=digits)
  }
  invisible(x)
}

as.data.frame.univaov <- function(x, row.names, optional, by=c("response", "term"), ...) {
  melt <- function(data, varnames = names(dimnames(data)), value.name = "value") {
    dn <- dimnames(data)
    labels <- expand.grid( dn[[1]], dn[[2]])
    colnames(labels) <- varnames
    value_df <- setNames(data.frame(as.vector(data)), value.name)
    cbind(labels, value_df)
  }
  nv <- ncol(x$F)
  nt <- nrow(x$F)
  by <- match.arg(by)
  if (by=="response") {
    vn <- c("term", "response")
    df <- matrix(x$SS[1:nt, "df", drop=FALSE], nrow=nt, ncol=nv)
    SS <- melt(x$SS[1:nt, -1, drop=FALSE], varnames=vn, value.name="SS")	
    F <- melt(x$F, varnames=vn, value.name="F")
    p <- melt(x$p, varnames=vn, value.name="p")
    if (!is.null(x$p.adjust)) p.adjust <- melt(x$p.adjust, varnames=vn, value.name="adjusted p")
  }
  else {
    vn <- rev(c("term", "response"))
    df <- t(matrix(x$SS[1:nt, "df", drop=FALSE], nrow=nt, ncol=nv))
    SS <- melt(t(x$SS[1:nt, -1, drop=FALSE]), varnames=vn, value.name="SS")	
    F <- melt(t(x$F), varnames=vn, value.name="F")
    p <- melt(t(x$p), varnames=vn, value.name="p")
    if (!is.null(x$p.adjust)) p.adjust <- melt(t(x$p.adjust), varnames=vn, value.name="adjusted p")
  }
  
  result <- cbind(SS[,c(2,1,3)], df=c(df), F=F[,"F"], p=p[,"p"])
  if (!is.null(x$p.adjust)) result <- cbind(result, "adjusted p"=p.adjust[, "adjusted p"])
  result
}

Anova.manova <- function(mod, ...){
  class(mod) <- c("mlm", "lm")
  Anova(mod, ...)
}

Manova <- function(mod, ...){
  UseMethod("Manova")
}

Manova.mlm <- function(mod, ...){
  Anova(mod, ...)
}

# Cox regression models

df.residual.coxph <- function(object, ...){    
  object$n - sum(!is.na(coef(object)))
}

alias.coxph <- function(object, ...){
  if(any(which <- is.na(coef(object)))) return(list(Complete=which))
  else list()
}

logLik.coxph <- function(object, ...) object$loglik[2]

Anova.coxph <- function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald"), ...){
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)
  if (length((mod$rscore) > 0) && (test.statistic == "LR")){ 
    warning("LR tests unavailable with robust variances\n  Wald tests substituted")
    test.statistic <- "Wald"
  }
  rhs <- as.character(formula(mod))[[3]]
  if (test.statistic == "LR" && grepl("tt\\(", rhs)){
    warning("LR tests unavailable for models using the tt argument\n  Wald tests substituted")
    test.statistic <- "Wald"
  }
  names <- term.names(mod)
  clusters <- grepl("cluster\\(", names)
  strata <- grepl("strata\\(", names)
  if ((any(clusters) || any(strata)) && test.statistic == "LR"){
    warning("LR tests not supported for models with clusters or strata\n Wald tests substituted")
    test.statistic <- "Wald"
  }
  switch(type,
         II=switch(test.statistic,
                   LR=Anova_II_LR_coxph(mod),
                   Wald=Anova.default(mod, type="II", test.statistic="Chisq", vcov.=vcov(mod, complete=FALSE))),
         III=switch(test.statistic,
                    LR=Anova_III_LR_coxph(mod),
                    Wald=Anova.default(mod, type="III", test.statistic="Chisq", vcov.=vcov(mod, complete=FALSE))),
         "2"=switch(test.statistic,
                    LR=Anova_II_LR_coxph(mod),
                    Wald=Anova.default(mod, type="II", test.statistic="Chisq", vcov.=vcov(mod, complete=FALSE))),
         "3"=switch(test.statistic,
                    LR=Anova_III_LR_coxph(mod),
                    Wald=Anova.default(mod, type="III", test.statistic="Chisq", vcov.=vcov(mod, complete=FALSE))))
}

Anova_II_LR_coxph <- function(mod, ...){
  if (!requireNamespace("survival")) stop("survival package is missing")
  which.nms <- function(name) which(asgn == which(names == name))
  fac <-attr(terms(mod), "factors")
  names <- term.names(mod)
  n.terms <- length(names)
  df <- df.terms(mod)
  if (sum(df > 0) < 2) {
    return(anova(mod, test="Chisq"))
  }
  method <- mod$method
  weights <- mod$weights
  X <- model.matrix(mod)
  asgn <- attr(X, 'assign')
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  for (term in 1:n.terms){
    if (df[names[term]] == 0){
      message("skipping term ", names[term])
      next
    }
    rels <- names[relatives(names[term], names, fac)]
    exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
    mod.1 <- survival::coxph(mod$y ~ X[, -exclude.1, drop = FALSE], method=method, weights=weights)
    loglik.1 <- logLik(mod.1)
    mod.2 <- if (length(rels) == 0) mod
    else {
      exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
      survival::coxph(mod$y ~ X[, -exclude.2, drop = FALSE], method=method, weights=weights)
    }
    loglik.2 <- logLik(mod.2)
    LR[term] <- -2*(loglik.1 - loglik.2)
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  result <- result[df > 0, , drop=FALSE]
  attr(result, "heading") <- "Analysis of Deviance Table (Type II tests)"
  result
}

Anova_III_LR_coxph <- function(mod, ...){
  if (!requireNamespace("survival")) stop("survival package is missing")
  which.nms <- function(name) which(asgn == which(names == name))
  fac <-attr(terms(mod), "factors")
  names <- term.names(mod)
  n.terms <- length(names)
  df <- df.terms(mod)
  if (sum(df > 0) < 2) {
    return(anova(mod, test="Chisq"))
  }
  method <- mod$method
  weights <- mod$weights
  X <- model.matrix(mod)
  asgn <- attr(X, 'assign')
  LR <- p <- rep(0, n.terms)
  loglik1 <- logLik(mod)
  for (term in 1:n.terms){
    if (df[names[term]] == 0){
      message("skipping term ", names[term])
      next
    }
    mod.0 <- survival::coxph(mod$y ~ X[, -which.nms(names[term])], method=method, weights=weights)
    LR[term] <- -2*(logLik(mod.0) - loglik1)
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df","Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  result <- result[df > 0, , drop=FALSE]
  attr(result,"heading") <- "Analysis of Deviance Table (Type III tests)"
  result
}

# parametric survival regression models

alias.survreg <- function(object, ...){
  if(any(which <- diag(vcov(object, complete=FALSE)) < 1e-10)) return(list(Complete=which))
  else list()
}

logLik.survreg <- function(object, ...) object$loglik[2]

Anova.survreg <- function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald"), ...){
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)
  if (length((mod$rscore) > 0) && (test.statistic == "LR")){ 
    warning("LR tests unavailable with robust variances\nWald tests substituted")
    test.statistic <- "Wald"
  }
  names <- term.names(mod)
  clusters <- grepl("cluster\\(", names)
  strata <- grepl("strata\\(", names)
  if ((any(clusters) || any(strata)) && test.statistic == "LR"){
    warning("LR tests not supported for models with clusters or strata\n Wald tests substituted")
    test.statistic <- "Wald"
  }
  switch(type,
         II=switch(test.statistic,
                   LR=Anova_II_LR_survreg(mod),
                   Wald=Anova_Wald_survreg(mod, type="2")),
         III=switch(test.statistic,
                    LR=Anova_III_LR_survreg(mod),
                    Wald=Anova_Wald_survreg(mod, type="3")),
         "2"=switch(test.statistic,
                    LR=Anova_II_LR_survreg(mod),
                    Wald=Anova_Wald_survreg(mod, type="2")),
         "3"=switch(test.statistic,
                    LR=Anova_III_LR_survreg(mod),
                    Wald=Anova_Wald_survreg(mod, type="3")))
}

Anova_II_LR_survreg <- function(mod, ...){
  if (!requireNamespace("survival")) stop("survival package is missing")
  dist <- mod$dist
  scale <- mod$call$scale
  weights <- model.frame(mod)$"(weights)"
  arg.list <- list(dist=dist)
  if (!is.null(scale)) arg.list$scale <- scale
  if (!is.null(weights)) arg.list$weights <- weights
  which.nms <- function(name) which(asgn == which(names == name))
  fac <-attr(terms(mod), "factors")
  names <- term.names(mod)
  X <- model.matrix(mod)
  asgn <- attr(X, 'assign')
  asgn <- asgn[asgn != 0]
  if (has_intercept(mod)){
    int <- which(names == "(Intercept)")
    X <- X[, -int]
    names <- names[-int]
  }
  n.terms <- length(names)
  if (n.terms < 2) {
    result <- anova(mod)
    if (!inherits(result, "anova")) class(result) <- c("anova", class(result))
    return(result)
  }
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  y <- model.frame(mod)[,1]
  for (term in 1:n.terms){
    rels <- names[relatives(names[term], names, fac)]
    exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
    arg.list$formula <- y ~ X[, -exclude.1, drop = FALSE]
    mod.1 <- do.call(survival::survreg, arg.list)
    loglik.1 <- logLik(mod.1)
    mod.2 <- if (length(rels) == 0) mod
    else {
      arg.list$formula <- y ~ X[, -exclude.2, drop = FALSE]
      exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
      do.call(survival::survreg, arg.list)     
    }
    loglik.2 <- logLik(mod.2)
    LR[term] <- -2*(loglik.1 - loglik.2)
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- "Analysis of Deviance Table (Type II tests)"
  result
}

Anova_III_LR_survreg <- function(mod, ...){
  if (!requireNamespace("survival")) stop("survival package is missing")
  dist <- mod$dist
  scale <- mod$call$scale
  weights <- model.frame(mod)$"(weights)"
  arg.list <- list(dist=dist)
  if (!is.null(scale)) arg.list$scale <- scale
  if (!is.null(weights)) arg.list$weights <- weights
  which.nms <- function(name) which(asgn == which(names == name))
  fac <-attr(terms(mod), "factors")
  names <- term.names(mod)
  X <- model.matrix(mod)
  asgn <- attr(X, 'assign')
  asgn <- asgn[asgn != 0]
  if (has_intercept(mod)){
    int <- which(names == "(Intercept)")
    X <- X[, -int]
    names <- names[-int]
  }
  n.terms <- length(names)
  if (n.terms < 2){
    result <- anova(mod)
    if (!inherits(result, "anova")) class(result) <- c("anova", class(result))
    return(result)
  }
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  y <- model.frame(mod)[,1]
  loglik1 <- logLik(mod)
  for (term in 1:n.terms){
    arg.list$formula <- y ~ X[, -which.nms(names[term])]
    mod.0 <- do.call(survival::survreg, arg.list)
    LR[term] <- -2*(logLik(mod.0) - loglik1)
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df","Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result,"heading") <- "Analysis of Deviance Table (Type III tests)"
  result
}

Anova_Wald_survreg <- function(mod, type){
  V <- vcov(mod, complete=FALSE)
  b <- coef(mod)
  if (length(b) != nrow(V)){
    p <- which(grepl("^Log\\(scale", rownames(V)))
    if (length(p) > 0) V <- V[-p, -p]
  }
  Anova.default(mod, V, test.statistic="Chisq", type=type)
}

# Default Anova() method: requires methods for vcov() (if vcov. argument not specified) and coef().

Anova.default <- function(mod, type=c("II","III", 2, 3), test.statistic=c("Chisq", "F"), 
                          vcov.=vcov(mod, complete=FALSE), singular.ok, error.df, ...){
  if (missing(error.df)){
    error.df <- df.residual(mod)
    test.statistic <- match.arg(test.statistic)
    if (test.statistic == "F" && (is.null(error.df) || is.na(error.df))){
      test.statistic <- "Chisq"
      message("residual df unavailable, test.statistic set to 'Chisq'")
    }
  }
  vcov. <- getVcov(vcov., mod)
  X <- model.matrix(mod)
  if (!is.matrix(X)) {
    warning("result of model.matrix(mod) is not a matrix",
            "\ntests may be incorrect\n")
  } else {
    coef.names <- colnames(X)
    if (any(bad <- !(names(coef(mod)) %in% coef.names)))
      warning("there are coefficients in coef(mod)",
              "\nthat are not in the model matrix:\n",
              paste(names(coef(mod))[bad], collapse=", "),
              "\ntests may be incorrect\n\n")
    if (any(bad <- !(colnames(vcov.) %in% coef.names))) 
      warning("there are rows/columns in vcov.",
              "\nthat are not in the model matrix:\n",
              paste(colnames(vcov.)[bad], collapse=", "),
              "\ntests may be incorrect")
  }
  type <- as.character(type)
  type <- match.arg(type)
  if (missing(singular.ok))
    singular.ok <- type == "2" || type == "II"
  switch(type,
         II=Anova_II_default(mod, vcov., test.statistic, singular.ok=singular.ok, 
                             error.df=error.df),
         III=Anova_III_default(mod, vcov., test.statistic, singular.ok=singular.ok,
                               error.df=error.df),
         "2"=Anova_II_default(mod, vcov., test.statistic, singular.ok=singular.ok,
                              error.df=error.df),
         "3"=Anova_III_default(mod, vcov., test.statistic, singular.ok=singular.ok,
                               error.df=error.df))
}

assignVector <- function(model, ...) UseMethod("assignVector")

assignVector.default <- function(model, ...){
  m <- model.matrix(model)
  assign <- attr(m, "assign")
  if (!is.null(assign)) {
    if (!has_intercept(model)) assign <- assign[assign != 0]
    return (assign)
  }
  m <- model.matrix(formula(model), data=model.frame(model))
  assign <- attr(m, "assign")
  if (!has_intercept(model)) assign <- assign[assign != 0]
  assign
}

Anova_II_default <- function(mod, vcov., test, singular.ok=TRUE, error.df,...){
  hyp.term <- function(term){
    which.term <- which(term==names)
    subs.term <- if (is.list(assign)) assign[[which.term]] else which(assign == which.term)
    relatives <- relatives(term, names, fac)
    subs.relatives <- NULL
    for (relative in relatives){
      sr <- if (is.list(assign)) assign[[relative]] else which(assign == relative)
      subs.relatives <- c(subs.relatives, sr)
    }
    hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
    hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
    hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
    hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
    hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
    else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov.))            
    hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
                                              function(x) all(x == 0)), , drop=FALSE]
    if (nrow(hyp.matrix.term) == 0)
      return(c(statistic=NA, df=0))            
    hyp <- linearHypothesis.default(mod, hyp.matrix.term, 
                                    vcov.=vcov., test=test, error.df=error.df,
                                    singular.ok=singular.ok, ...)
    if (test=="Chisq") c(statistic=hyp$Chisq[2], df=hyp$Df[2])
    else c(statistic=hyp$F[2], df=hyp$Df[2])
  }
  not.aliased <- !is.na(coef(mod))
  if (!singular.ok && !all(not.aliased))
    stop("there are aliased coefficients in the model")
  fac <- attr(terms(mod), "factors")
  intercept <- has_intercept(mod)
  p <- length(coefficients(mod))
  I.p <- diag(p)
  assign <- assignVector(mod) 
  if (!is.list(assign)) assign[!not.aliased] <- NA
  else if (intercept) assign <- assign[-1]
  names <- term.names(mod)
  if (intercept) names <- names[-1]
  n.terms <- length(names)
  df <- c(rep(0, n.terms), error.df)
  if (inherits(mod, "coxph") || inherits(mod, "survreg")){
    if (inherits(mod, "coxph")) assign <- assign[assign != 0]
    clusters <- grep("^cluster\\(", names)
    strata <- grep("^strata\\(.*\\)$", names)
    for (cl in clusters) assign[assign > cl] <- assign[assign > cl] - 1
    for (st in strata) assign[assign > st] <- assign[assign > st] - 1
    if (length(clusters) > 0 || length(strata) > 0) {
      message("skipping term ", paste(names[c(clusters, strata)], collapse=", "))
      names <- names[-c(clusters, strata)]
      df <- df[-c(clusters, strata)]
      n.terms <- n.terms - length(clusters) - length(strata)
    }
  }
  p <- teststat <- rep(0, n.terms + 1)
  teststat[n.terms + 1] <- p[n.terms + 1] <- NA
  for (i in 1:n.terms){
    hyp <- hyp.term(names[i])
    teststat[i] <- abs(hyp["statistic"])
    df[i] <- abs(hyp["df"])
    p[i] <- if (test == "Chisq") 
      pchisq(teststat[i], df[i], lower.tail=FALSE) 
    else pf(teststat[i], df[i], df[n.terms + 1], lower.tail=FALSE)
  }
  result <- if (test == "Chisq"){ 
    if (length(df) == n.terms + 1) df <- df[1:n.terms]
    data.frame(df[df > 0], teststat[!is.na(teststat)], p[!is.na(teststat)])
  }
  else data.frame(df, teststat, p)
  if (nrow(result) == length(names) + 1) names <- c(names,"Residuals")
  row.names(result) <- names[df > 0]
  names(result) <- c ("Df", test, if (test == "Chisq") "Pr(>Chisq)" 
                      else "Pr(>F)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
                               paste("Response:", responseName(mod)))
  result
}

Anova_III_default <- function(mod, vcov., test, singular.ok=FALSE, error.df, ...){
  intercept <- has_intercept(mod)
  p <- length(coefficients(mod))
  I.p <- diag(p)
  names <- term.names(mod)
  n.terms <- length(names)
  assign <- assignVector(mod) 
  df <- c(rep(0, n.terms), error.df)
  if (inherits(mod, "coxph")){
    if (intercept) names <- names[-1]
    assign <- assign[assign != 0]
    clusters <- grep("^cluster\\(", names)
    strata <- grep("^strata\\(.*\\)$", names)
    for (cl in clusters) assign[assign > cl] <- assign[assign > cl] - 1
    for (st in strata) assign[assign > st] <- assign[assign > st] - 1
    if (length(clusters) > 0 || length(strata) > 0) {
      message("skipping term ", paste(names[c(clusters, strata)], collapse=", "))
      names <- names[-c(clusters, strata)]
      df <- df[-c(clusters, strata)]
      n.terms <- n.terms - length(clusters) - length(strata)
    }
  }
  if (intercept) df[1] <- sum(grepl("^\\(Intercept\\)", names(coef(mod))))
  teststat <- rep(0, n.terms + 1)
  p <- rep(0, n.terms + 1)
  teststat[n.terms + 1] <- p[n.terms + 1] <- NA
  not.aliased <- !is.na(coef(mod))
  if (!singular.ok && !all(not.aliased))
    stop("there are aliased coefficients in the model")
  for (term in 1:n.terms){
    subs <- if (is.list(assign)) assign[[term]] else which(assign == term - intercept)    
    hyp.matrix <- I.p[subs,,drop=FALSE]
    hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
    hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
    if (nrow(hyp.matrix) == 0){
      teststat[term] <- NA
      df[term] <- 0
      p[term] <- NA
    }
    else {
      hyp <- linearHypothesis.default(mod, hyp.matrix, 
                                      vcov.=vcov., test=test, error.df=error.df,
                                      singular.ok=singular.ok, ...)
      teststat[term] <- if (test=="Chisq") hyp$Chisq[2] else hyp$F[2]
      df[term] <- abs(hyp$Df[2])
      p[term] <- if (test == "Chisq") 
        pchisq(teststat[term], df[term], lower.tail=FALSE) 
      else pf(teststat[term], df[term], df[n.terms + 1], lower.tail=FALSE)
    }
  }
  result <- if (test == "Chisq"){ 
    if (length(df) == n.terms + 1) df <- df[1:n.terms]
    data.frame(df, teststat[!is.na(teststat)], p[!is.na(teststat)])
  }
  else data.frame(df, teststat, p)
  if (nrow(result) == length(names) + 1) names <- c(names,"Residuals")
  row.names(result) <- names
  names(result) <- c ("Df", test, if (test == "Chisq") "Pr(>Chisq)" 
                      else "Pr(>F)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", 
                               paste("Response:", responseName(mod)))
  result
}

## functions for mixed models

# the following function, not exported, to make car consistent with CRAN and development versions of lme4 and with nlme

fixef <- function (object){
  if (isS4(object)) {
    if (!inherits(object, "merMod")) object@fixef
    else lme4::fixef(object)
  }
  else object$coefficients$fixed
}

Anova.merMod <- function(mod, type=c("II","III", 2, 3), 
                         test.statistic=c("Chisq", "F"),
                         vcov.=vcov(mod, complete=FALSE), singular.ok, ...){
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)
  if (!missing(vcov.)) {
    if (test.statistic != "F"){
      message("Coefficient covariances computed by ", deparse(substitute(vcov.)))
    } else {
      warning('test.statistic="F"; vcov. argument ignored')
    }
  }
  vcov. <- getVcov(vcov., mod)
  if (missing(singular.ok))
    singular.ok <- type == "2" || type == "II"
  Anova.mer(mod=mod, type=type, test.statistic=test.statistic, vcov.=vcov.,
            singular.ok=singular.ok, ...)
}

Anova.mer <- function(mod, type=c("II","III", 2, 3), test.statistic=c("Chisq", "F"),
                      vcov.=vcov(mod, complete=FALSE), singular.ok, ...){
  vcov. <- getVcov(vcov., mod)
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)
  if (missing(singular.ok))
    singular.ok <- type == "2" || type == "II"
  switch(type,
         II=Anova_II_mer(mod, test=test.statistic, vcov., singular.ok=singular.ok),
         III=Anova_III_mer(mod, test=test.statistic,  vcov., singular.ok=singular.ok),
         "2"=Anova_II_mer(mod, test=test.statistic, vcov., singular.ok=singular.ok),
         "3"=Anova_III_mer(mod, test=test.statistic, vcov., singular.ok=singular.ok))
}

Anova_II_mer <- function(mod, vcov., singular.ok=TRUE, test=c("Chisq", "F"), ...){
  hyp.term <- function(term){ 
    which.term <- which(term==names)
    subs.term <- which(assign==which.term)
    relatives <- relatives(term, names, fac)
    subs.relatives <- NULL
    for (relative in relatives) 
      subs.relatives <- c(subs.relatives, which(assign==relative))
    hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
    hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
    hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
    hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
    hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
    else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov.))            
    hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
                                              function(x) all(x == 0)), , drop=FALSE]
    if (nrow(hyp.matrix.term) == 0)
      return(c(statistic=NA, df=0))            
    hyp <- linearHypothesis(mod, hyp.matrix.term, 
                            vcov.=vcov., singular.ok=singular.ok, test=test, ...)
    if (test == "Chisq") return(c(statistic=hyp$Chisq[2], df=hyp$Df[2]))
    else return(c(statistic=hyp$F[2], df=hyp$Df[2], res.df=hyp$Res.Df[2]))
  }
  test <- match.arg(test)
  not.aliased <- !is.na(fixef(mod))
  if (!singular.ok && !all(not.aliased))
    stop("there are aliased coefficients in the model")
  fac <- attr(terms(mod), "factors")
  intercept <- has_intercept(mod)
  p <- length(fixef(mod))
  I.p <- diag(p)
  if (test == "F"){
    vcov. <- as.matrix(pbkrtest::vcovAdj(mod, details=0))
  }
  assign <- attr(model.matrix(mod), "assign")
  assign[!not.aliased] <- NA
  names <- term.names(mod)
  if (intercept) names <- names[-1]
  n.terms <- length(names)
  p <- teststat <- df <- res.df <- rep(0, n.terms)
  for (i in 1:n.terms){
    hyp <- hyp.term(names[i])
    teststat[i] <- abs(hyp["statistic"])
    df[i] <- abs(hyp["df"])
    res.df[i] <- hyp["res.df"]
    p[i] <- if (test == "Chisq") pchisq(teststat[i], df[i], lower.tail=FALSE) 
    else pf(teststat[i], df[i], res.df[i], lower.tail=FALSE)
  } 
  if (test=="Chisq"){
    result <- data.frame(teststat, df, p)
    row.names(result) <- names
    names(result) <- c ("Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type II Wald chisquare tests)\n", 
                                 paste("Response:", responseName(mod)))
  }
  else {
    result <- data.frame(teststat, df, res.df, p)
    row.names(result) <- names
    names(result) <- c ("F", "Df", "Df.res", "Pr(>F)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type II Wald F tests with Kenward-Roger df)\n", 
                                 paste("Response:", responseName(mod)))
  }
  result
}

Anova_III_mer <- function(mod, vcov., singular.ok=FALSE, test=c("Chisq", "F"), ...){
  intercept <- has_intercept(mod)
  p <- length(fixef(mod))
  I.p <- diag(p)
  names <- term.names(mod)
  n.terms <- length(names)
  assign <- attr(model.matrix(mod), "assign")
  p <- teststat <- df <- res.df <- rep(0, n.terms)
  if (intercept) df[1] <- 1
  not.aliased <- !is.na(fixef(mod))
  if (!singular.ok && !all(not.aliased))
    stop("there are aliased coefficients in the model")
  if (test == "F"){
    vcov. <- as.matrix(pbkrtest::vcovAdj(mod, details=0))
  }
  for (term in 1:n.terms){
    subs <- which(assign == term - intercept)        
    hyp.matrix <- I.p[subs,,drop=FALSE]
    hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
    hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
    if (nrow(hyp.matrix) == 0){
      teststat[term] <- NA
      df[term] <- 0
      p[term] <- NA
    }
    else {
      hyp <- linearHypothesis(mod, hyp.matrix, test=test,
                              vcov.=vcov., singular.ok=singular.ok, ...)
      if (test == "Chisq"){
        teststat[term] <-  hyp$Chisq[2] 
        df[term] <- abs(hyp$Df[2])
        p[term] <- pchisq(teststat[term], df[term], lower.tail=FALSE) 
      }
      else{
        teststat[term] <-  hyp$F[2]
        df[term] <- abs(hyp$Df[2])
        res.df[term]=hyp$Res.Df[2]
        p[term] <- pf(teststat[term], df[term], res.df[term], lower.tail=FALSE)
      }
    }
  }
  if (test == "Chisq"){
    result <- data.frame(teststat, df, p)
    row.names(result) <- names
    names(result) <- c ("Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type III Wald chisquare tests)\n", 
                                 paste("Response:", responseName(mod)))
  }
  else {
    result <- data.frame(teststat, df, res.df, p)
    row.names(result) <- names
    names(result) <- c ("F", "Df", "Df.res", "Pr(>F)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Analysis of Deviance Table (Type III Wald F tests with Kenward-Roger df)\n", 
                                 paste("Response:", responseName(mod)))
  }
  result
}

has_intercept.lme <- function(model, ...){
  any(names(fixef(model)) == "(Intercept)")
}

Anova.lme <- function(mod, type=c("II","III", 2, 3),
                      vcov.=vcov(mod, complete=FALSE), singular.ok, ...){
  if (!missing(vcov.)) message("Coefficient covariances computed by ", deparse(substitute(vcov.)))
  vcov. <- getVcov(vcov., mod)
  type <- as.character(type)
  type <- match.arg(type)
  if (missing(singular.ok))
    singular.ok <- type == "2" || type == "II"
  switch(type,
         II=Anova_II_lme(mod, vcov., singular.ok=singular.ok),
         III=Anova_III_lme(mod, vcov., singular.ok=singular.ok),
         "2"=Anova_II_lme(mod, vcov., singular.ok=singular.ok),
         "3"=Anova_III_lme(mod, vcov., singular.ok=singular.ok))
}

Anova_II_lme <- function(mod, vcov., singular.ok=TRUE, ...){
  hyp.term <- function(term){
    which.term <- which(term==names)
    subs.term <- which(assign==which.term)
    relatives <- relatives(term, names, fac)
    subs.relatives <- NULL
    for (relative in relatives) 
      subs.relatives <- c(subs.relatives, which(assign==relative))
    hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
    hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
    hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
    hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]       
    hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
    else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov.))            
    hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
                                              function(x) all(x == 0)), , drop=FALSE]
    if (nrow(hyp.matrix.term) == 0)
      return(c(statistic=NA, df=0))            
    hyp <- linearHypothesis(mod, hyp.matrix.term, 
                            vcov.=vcov., singular.ok=singular.ok, ...)
    c(statistic=hyp$Chisq[2], df=hyp$Df[2])
  }
  not.aliased <- !is.na(fixef(mod))
  if (!singular.ok && !all(not.aliased))
    stop("there are aliased coefficients in the model")
  fac <- attr(terms(mod), "factors")
  intercept <- has_intercept(mod)
  p <- length(fixef(mod))
  I.p <- diag(p)
  attribs.mm <- attributes(model.matrix(mod))
  assign <- attribs.mm$assign
  nms.coef <- names(coef(mod))
  nms.mm <- attribs.mm$dimnames[[2]]
  assign[!not.aliased] <- NA
  valid.coefs <- nms.mm %in% nms.coef
  if (any(!valid.coefs)){
    warning("The following coefficients are not in the model due to missing levels:\n",
            paste(nms.mm[!valid.coefs], collapse=", "))
  }
  assign <- assign[valid.coefs]
  names <- term.names(mod)
  if (intercept) names <- names[-1]
  n.terms <- length(names)
  p <- teststat <- df <- rep(0, n.terms)
  for (i in 1:n.terms){
    hyp <- hyp.term(names[i])
    teststat[i] <- abs(hyp["statistic"])
    df[i] <- abs(hyp["df"])
    p[i] <- pchisq(teststat[i], df[i], lower.tail=FALSE) 
  }    
  result <- data.frame(teststat, df, p)
  row.names(result) <- names
  names(result) <- c("Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type II tests)\n", 
                               paste("Response:", responseName(mod)))
  result
}

Anova_III_lme <- function(mod, vcov., singular.ok=FALSE, ...){
  intercept <- has_intercept(mod)
  p <- length(fixef(mod))
  I.p <- diag(p)
  names <- term.names(mod)
  n.terms <- length(names)
  attribs.mm <- attributes(model.matrix(mod))
  assign <- attribs.mm$assign
  nms.coef <- names(coef(mod))
  nms.mm <- attribs.mm$dimnames[[2]]
  not.aliased <- !is.na(fixef(mod))
  if (!singular.ok && !all(not.aliased))
    stop("there are aliased coefficients in the model")
  assign[!not.aliased] <- NA
  valid.coefs <- nms.mm %in% nms.coef
  if (any(!valid.coefs)){
    warning("The following coefficients are not in the model due to missing levels:\n",
            paste(nms.mm[!valid.coefs], collapse=", "))
  }
  assign <- assign[valid.coefs]
  df <- rep(0, n.terms)
  if (intercept) df[1] <- 1
  p <- teststat <-rep(0, n.terms)
  for (term in 1:n.terms){
    subs <- which(assign == term - intercept)        
    hyp.matrix <- I.p[subs,,drop=FALSE]
    hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
    hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]        
    if (nrow(hyp.matrix) == 0){
      teststat[term] <- NA
      df[term] <- 0
      p[term] <- NA
    }
    else {
      hyp <- linearHypothesis(mod, hyp.matrix, 
                              vcov.=vcov., singular.ok=singular.ok, ...)
      teststat[term] <-  hyp$Chisq[2] 
      df[term] <- abs(hyp$Df[2])
      p[term] <- pchisq(teststat[term], df[term], lower.tail=FALSE) 
    }
  }
  result <- data.frame(teststat, df, p)
  row.names(result) <- names
  names(result) <- c ("Chisq",  "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- c("Analysis of Deviance Table (Type III tests)\n", 
                               paste("Response:", responseName(mod)))
  result
}

Anova.svyglm <- function(mod, ...) Anova.default(mod, ...)

Anova.rlm <- function(mod, ...) Anova.default(mod, test.statistic="F", ...)

Anova.coxme <- function(mod, type=c("II","III", 2, 3), test.statistic=c("Wald", "LR"), ...){
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)
  switch(type,
         II=switch(test.statistic,
                   LR=Anova_II_LR_coxme(mod, ...),
                   Wald=Anova.default(mod, type="II", test.statistic="Chisq", ...)),
         III=switch(test.statistic,
                    LR=stop("type-III LR tests not available for coxme models"),
                    Wald=Anova.default(mod, type="III", test.statistic="Chisq", ...)),
         "2"=switch(test.statistic,
                    LR=Anova_II_LR_coxme(mod, ...),
                    Wald=Anova.default(mod, type="II", test.statistic="Chisq", ...)),
         "3"=switch(test.statistic,
                    LR=stop("type-III LR tests not available for coxme models"),
                    Wald=Anova.default(mod, type="III", test.statistic="Chisq")))
}

Anova_II_LR_coxme <- function(mod, ...){
  if (!requireNamespace("coxme")) stop("coxme package is missing")
  which.nms <- function(name) which(asgn == which(names == name))
  fac <-attr(terms(mod), "factors")
  names <- term.names(mod)
  n.terms <- length(names)
  if (n.terms < 2){
    return(anova(mod, test="Chisq"))
  }
  X <- model.matrix(mod)
  asgn <- attr(X, 'assign')
  p <- LR <- rep(0, n.terms)
  df <- df.terms(mod)
  random <- mod$formulaList$random
  random <- sapply(random, as.character)[2, ]
  random <- paste(paste0("(", random, ")"), collapse=" + ")
  fixed <- as.character(mod$formulaList$fixed)[3]
  for (term in 1:n.terms){
    rels <- names[relatives(names[term], names, fac)]
    formula <- paste0(". ~ . - ", paste(c(names[term], rels), collapse=" - "), " + ", random)
    mod.1 <- update(mod, as.formula(formula))
    loglik.1 <- logLik(mod.1, type="integrated")
    mod.2 <- if (length(rels) == 0) mod
    else {
      formula <- paste0(". ~ . - ", paste(rels, collapse=" - "), " + ", random)
      update(mod, as.formula(formula))
    }
    loglik.2 <- logLik(mod.2, type="integrated")
    LR[term] <- -2*(loglik.1 - loglik.2)
    p[term] <- pchisq(LR[term], df[term], lower.tail=FALSE)
  }
  result <- data.frame(LR, df, p)
  row.names(result) <- names
  names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
  class(result) <- c("anova", "data.frame")
  attr(result, "heading") <- "Analysis of Deviance Table (Type II tests)"
  result
}

# the following unexported methods make Anova.default() and linearHypotheis.default() work with "svyolr" objects

assignVector.svyolr <- function(model, ...){
  m <- model.matrix(model)
  assign <- attr(m, "assign")
  assign[assign != 0]
}

coef.svyolr <- function(object, ...) NextMethod()

vcov.svyolr <- function(object, ...){
  nms <- names(coef(object))
  (object$var)[nms, nms]
}

# Anova() method for svycoxph objects, to prevent test="LR"

Anova.svycoxph <- function(mod, type=c("II", "III", 2, 3),
                           test.statistic="Wald", ...){
  test.statistic <- match.arg(test.statistic)
  type <- as.character(type)
  type <- match.arg(type)
  NextMethod(test.statistic=test.statistic, type=type, ...)
}

# Anova() methods for "clm" and "clmm" objects (ordinal package)
#   coef(), vcov(), and model.matrix() methods not exported

Anova.clm <- function(mod, ...){
  class(mod) <- c("clmAnova", class(mod))
  Anova.default(mod, ...)
}

coef.clmAnova <- function(object, ...) object$beta

vcov.clmAnova <- function(object, ...){
  coef.names <- names(coef(object))
  V <- NextMethod()
  V[coef.names, coef.names]
}

model.matrix.clmAnova <- function(object, ...){
  X <- NextMethod()
  X$X
}

Anova.clmm <- function(mod, ...){
  class(mod) <- c("clmmAnova", class(mod))
  Anova.default(mod, ...)
}


coef.clmmAnova <- function(object, ...) {
  names.thresholds <- colnames(object$Theta)
  coefs <- NextMethod()
  names.all <- names(coefs)
  names.fixed <- names.all[!names.all %in% names.thresholds]
  coefs[names.fixed]
}

vcov.clmmAnova <- function(object, ...){
  coef.names <- names(coef(object))
  V <- NextMethod() 
  V[coef.names, coef.names]
}
