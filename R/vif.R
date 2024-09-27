#-------------------------------------------------------------------------------
# Revision history:
# 2009-09-28 by J. Fox 
# 2013-05-21 replaced vif.lm with vif.default and added
#            model.matrix.gls to make gls models work. J. Fox
# 2015-01-13: fixed model.matrix.gls to work with models with formulas as object. J. Fox
# 2020-12-19: new polr and svyolr methods for ordinal regression models. J. Fox
# 2022-03-11: added new vif.lm() methods that handles interactions. J. Fox
# 2022-06-07: rework vif.lm() for interations. J. Fox
# 2024-05-14: has.intercept() -> has_intercept(). J. Fox
#-------------------------------------------------------------------------------

# Generalized Variance-Inflation Factors (John Fox and Henric Nilsson)

vif<-function(mod, ...){
	UseMethod("vif")
}

vif.default <- function(mod, ...) {
    if (any(is.na(coef(mod)))) 
        stop ("there are aliased coefficients in the model")
    v <- vcov(mod)
    assign <- attr(model.matrix(mod), "assign")
    if (names(coefficients(mod)[1]) == "(Intercept)") {
        v <- v[-1, -1]
        assign <- assign[-1]
    }
    else warning("No intercept: vifs may not be sensible.")
    terms <- labels(terms(mod))
    n.terms <- length(terms)
    if (n.terms < 2) stop("model contains fewer than 2 terms")
    R <- cov2cor(v)
    detR <- det(R)
    result <- matrix(0, n.terms, 3)
    rownames(result) <- terms
    colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
    for (term in 1:n.terms) {
        subs <- which(assign == term)
        result[term, 1] <- det(as.matrix(R[subs, subs])) *
            det(as.matrix(R[-subs, -subs])) / detR
        result[term, 2] <- length(subs)
    }
    if (all(result[, 2] == 1)) result <- result[, 1]
    else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
    result
}

vif.merMod <- function(mod, ...) {
  if (any(is.na(fixef(mod)))) 
    stop ("there are aliased coefficients in the model")
  v <- as.matrix(vcov(mod))
  assign <- attr(model.matrix(mod), "assign")
  if (names(fixef(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  }
  else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("model contains fewer than 2 terms")
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) *
      det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) result <- result[, 1]
  else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  result
}

model.matrix.gls <- function(object, ...){
    model.matrix(formula(object), data=eval(object$call$data))
}

vif.polr <- function(mod, ...) {
  if (any(is.na(coef(mod)))) 
    stop ("there are aliased coefficients in the model")
  v <- vcov(mod)
  nms <- names(coef(mod))
  v <- v[nms, nms]
  assign <- attr(model.matrix(mod), "assign")
  assign <- assign[assign != 0]
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("model contains fewer than 2 terms")
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) *
      det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) result <- result[, 1]
  else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  result
}

vif.svyolr <- function(mod, ...) {
  if (any(is.na(coef(mod)))) 
    stop ("there are aliased coefficients in the model")
  v <- vcov(mod)
  nms <- names(coef(mod, intercepts=FALSE))
  v <- v[nms, nms]
  assign <- attr(model.matrix(mod), "assign")
  assign <- assign[assign != 0]
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("model contains fewer than 2 terms")
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) *
      det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) result <- result[, 1]
  else result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  result
}

vif.lm <- function(mod, type=c("terms", "predictor"), ...){ 
  
  type <- match.arg(type)
  
  if (any(attr(terms(mod), "order") > 1)){
    if (type == "terms"){
      message("there are higher-order terms (interactions) in this model\n",
              "consider setting type = 'predictor'; see ?vif")
    }
  }
  
  weighted <- any(weights(mod) != 1)
  if ((inherits(mod, "glm") || weighted) && type != "terms"){
    warning("type = 'predictor' is available only for unweighted linear models;\n",
            "  type = 'terms' will be used")
  }
  if (type == "terms" || weighted || inherits(mod, "glm")) {
    return(NextMethod())
  }
  
  factors <- attr(terms(mod), "factors")
  names <- term.names(mod)
  X <- model.matrix(mod)
  intercept <- has_intercept(mod)
  if (intercept) {
    names <- names[-1]
    assign.X <- attr(X, "assign")[-1]
    X <- X[, -1]
  } else {
    warning("No intercept: (G)VIFs may not be sensible.")
  }
  R <- cor(X)
  detR <- det(R)
  X.names <- colnames(X)
  formula <- formula(mod)[-2]
  terms <- attr(terms(mod), "term.labels")
  term.vars <- lapply(parse(text=terms), all.vars)
  predictors <- all.vars(formula)
  vifs <- matrix(0, length(predictors), 3)
  rownames(vifs) <- predictors
  colnames(vifs) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
  vifs <- as.data.frame(vifs)
  vifs$`Interacts With` <- rep("", length(predictors))
  vifs$`Other Predictors` <- rep("", length(predictors))
  all.cols <- 1:ncol(X)
  for (predictor in predictors){
    which.terms <- sapply(term.vars, function(vars) predictor %in% vars)
    related <- unique(unlist(strsplit(paste(terms[which.terms], collapse=":"), ":")))
    vifs[predictor, 4] <- if (length(related[-1]) > 0) {
      paste(related[-1], collapse=", ")
    } else {
      "--  "
    }
    unrelated <- setdiff(predictors, related)
    if (length(unrelated) > 0){
      unrelated.terms <- sapply(term.vars, 
                                function(vars) unrelated %in% vars)
      if (is.matrix(unrelated.terms)) unrelated.terms <- apply(unrelated.terms, 2, any)
      columns <- setdiff(all.cols, which(assign.X %in% which(unrelated.terms)))
      gvif <- det(R[columns, columns, drop=FALSE])*det(R[-columns, -columns, drop=FALSE])/detR
      vifs[predictor, 5] <- paste(unrelated, collapse=", ")
    } else {
      columns <- all.cols
      gvif <- 1
      vifs[predictor, 5] <- "--  "
    }
    p <- length(columns)
    vifs[predictor, 1:3] <- c(gvif, p, gvif^(1/(2*p)))
  }
  
  if (all(vifs[, 2] == 1)) {
    message("VIFs computed for predictors")
    return(vifs[, 1])
  } else {
    message("GVIFs computed for predictors")
    return(vifs)
  }
}
