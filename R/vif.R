#-------------------------------------------------------------------------------
# Revision history:
# 2009-09-28 by J. Fox 
# 2013-05-21 replaced vif.lm with vif.default and added
#            model.matrix.gls to make gls models work. J. Fox
# 2015-01-13: fixed model.matrix.gls to work with models with formulas as object. J. Fox
# 2020-12-19: new polr and svyolr methods for ordinal regression models. J. Fox
# 2022-03-11: added new vif.lm() methods that handles interactions. J. Fox
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

vif.lm <- function(mod, type=c("terms", "marginal", "high-order"), ...){ 
  
  highOrderTerms <- function(){
    high <- character(0)
    for (term in names){
      if (0 == length(relatives(term, names, factors))) high <- c(high, term)
    }
    high
  }
  
  lowerOrderRelatives <- function(term){
    candidates <- names[apply(factors, 2, function(x) any(factors[, term] & x))]
    first.order <- names[colSums(factors) == 1]
    if (term %in% first.order) return(NULL)
    excluded <- first.order[
      !sapply(first.order, function(x) term %in% names[relatives(x, names, factors)])
    ]
    lower <- if (length(excluded) == 0) candidates else candidates[factors[excluded, candidates] == 0]
    list(candidates=candidates, first.order=first.order, excluded=excluded, lower=lower)
    lower[lower != term]
  }
  
  toCompareHighorder <- function(term){
    relatives <- lowerOrderRelatives(term)
    which.relatives <- assign.X %in% sapply(c(term, relatives), function(x)  which(x == names))
    x1 <- X.names[which.relatives]
    x2 <- setdiff(X.names, x1)
    list(x1=x1, x2=x2)
  }
  
  toCompare <- function(term){
    relatives <- relatives(term, names, factors)
    x1.x2 <- X.names[!(assign.X %in% relatives)]
    x1 <- X.names[which(term == names) == assign.X]
    x2 <- setdiff(x1.x2, x1)
    list(x1=x1, x2=x2, x1.x2=x1.x2)
  }
  
  type <- match.arg(type)
  
  if (any(attr(terms(mod), "order") > 1)){
    if (type == "terms"){
      message("there are higher-order terms (interactions) in this model\n",
              "consider setting terms = 'marginal' or 'high-order'; see ?vif")
    }
  }
  
  weighted <- any(weights(mod) != 1)
  if ((inherits(mod, "glm") || weighted) && type != "terms"){
    warning("type = '", type, "' is available only for unweighted linear models;\n",
            "  type = 'terms' will be used")
  }
  if (type == "terms" || weighted || inherits(mod, "glm")) {
    return(NextMethod())
  }
  
  factors <- attr(terms(mod), "factors")
  names <- term.names(mod)
  X <- model.matrix(mod)
  intercept <- has.intercept(mod)
  if (intercept) {
    names <- names[-1]
    assign.X <- attr(X, "assign")[-1]
    X <- X[, -1]
  } else {
    warning("No intercept: (G)VIFs may not be sensible.")
  }
  R <- cor(X)
  X.names <- colnames(X)
  
  if (type == "marginal"){
    n.terms <- length(names)
    vifs <- matrix(0, n.terms, 3)
    colnames(vifs) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
    rownames(vifs) <- names
    for (i in 1:n.terms){
      compare <- toCompare(names[i])
      vifs[i, 1] <- det(R[compare$x1, compare$x1, drop=FALSE])*
        det(R[compare$x2, compare$x2, drop=FALSE])/
        det(R[compare$x1.x2, compare$x1.x2, drop=FALSE])
      p <- length(compare$x1)
      vifs[i, 2] <- p
      vifs[i, 3] <- vifs[i, 1]^(1/(2*p))
    }
    
  } else {
    high.order.terms <- highOrderTerms()
    n.terms <- length(high.order.terms)
    vifs <- matrix(0, n.terms, 3)
    colnames(vifs) <- c("GVIF", "Df", "GVIF^(1/(2*Df))")
    detR <- det(R)
    for (i in 1:n.terms){
      compare <- toCompareHighorder(high.order.terms[i])
      vifs[i, 1] <- det(R[compare$x1, compare$x1, drop=FALSE])*
        det(R[compare$x2, compare$x2, drop=FALSE])/
        detR
      p <- length(compare$x1)
      vifs[i, 2] <- p
      vifs[i, 3] <- vifs[i, 1]^(1/(2*p))
    }
    rownames(vifs) <- gsub("\\:", "\\*", high.order.terms)
  }
  
  if (all(vifs[, 2] == 1)) {
    if (type == "marginal") message("VIFs computed respecting marginality")
    else message("VIFs computed for high-order terms")
    return(vifs[, 1])
  } else {
    if (type == "marginal") message("GVIFs computed respecting marginality")
    else message("GVIFs computed for high-order terms")
    return(vifs)
  }
}
