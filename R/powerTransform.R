# 2009-09-16: added ... argument to print.summary.powerTransform. J. Fox
# 2015-02-02: added 'gamma' argument to get transformation of (U + gamma)
# 2015-08-10: added estimateTransform as a generic function
# 2015-08-24: made 'family' an explicit argument to powerTransformation to clairfy man page.
# 2017-01-28: bug-fix in yjPower
# 2017-05-02: function updates to accomodate bcnPower family.  S. Weisberg
# 2017-05-19: Changed summary.powerTransform; deleted invalid test; added roundlam to output
# 2017-07-17:  Added family object in return of estimateTransform.default; changed print function of summary.powerTransform B. Price
# 2017-10-25: modified print.powerTransform() and print.summary.powerTransform()
#             so that singular words are used for 1 parameter (e.g., "is" vs "are"). J. Fox
# 2017-12-01: removed plot.powerTransform


### Power families:
basicPower <- function(U,lambda, gamma=NULL) {
 if(!is.null(gamma)) basicPower(t(t(as.matrix(U) + gamma)), lambda) else{
 bp1 <- function(U,lambda){
  if(any(U[!is.na(U)] <= 0)) stop("First argument must be strictly positive.")
  if (abs(lambda) <= 1.e-6) log(U) else (U^lambda)
  }
  out <- U
  out <- if(is.matrix(out) | is.data.frame(out)){
    if(is.null(colnames(out))) colnames(out) <-
        paste("Z", 1:dim(out)[2],sep="")
    for (j in 1:ncol(out)) {out[, j] <- bp1(out[, j],lambda[j])
       colnames(out)[j] <- if(abs(lambda[j]) <= 1.e-6)
           paste("log(", colnames(out)[j],")", sep="") else
           paste(colnames(out)[j], round(lambda[j], 2), sep="^")}
    out}  else
    bp1(out, lambda)
  out}}

bcPower <- function(U, lambda, jacobian.adjusted=FALSE, gamma=NULL) {
 if(!is.null(gamma)) bcPower(t(t(as.matrix(U) + gamma)), lambda, jacobian.adjusted) else{
 bc1 <- function(U, lambda){
  if(any(U[!is.na(U)] <= 0)) stop("First argument must be strictly positive.")
  z <- if (abs(lambda) <= 1.e-6) log(U) else ((U^lambda) - 1)/lambda
  if (jacobian.adjusted == TRUE) {
    z * (exp(mean(log(U), na.rm=TRUE)))^(1-lambda)} else z
  }
  out <- U
  out <- if(is.matrix(out) | is.data.frame(out)){
    if(is.null(colnames(out))) colnames(out) <-
        paste("Z", 1:dim(out)[2], sep="")
    for (j in 1:ncol(out)) {out[, j] <- bc1(out[, j], lambda[j]) }
    colnames(out) <- paste(colnames(out), round(lambda, 2), sep="^")
    out}  else
    bc1(out, lambda)
  out}}

yjPower <- function(U, lambda, jacobian.adjusted=FALSE) {
 yj1 <- function(U, lambda){
  nonnegs <- U >= 0
  z <- rep(NA, length(U))
  z[which(nonnegs)] <- bcPower(U[which(nonnegs)]+1, lambda, jacobian.adjusted=FALSE)
  z[which(!nonnegs)] <- -bcPower(-U[which(!nonnegs)]+1, 2-lambda, jacobian.adjusted=FALSE)
  if (jacobian.adjusted == TRUE)
        z * (exp(mean(log((1 + abs(U))^(2 * nonnegs - 1)), na.rm=TRUE)))^(1 -
            lambda)
    else z
  }
  out <- U
  out <- if(is.matrix(out) | is.data.frame(out)){
    if(is.null(colnames(out))) colnames(out) <-
        paste("Z", 1:dim(out)[2], sep="")
    for (j in 1:ncol(out)) {out[, j] <- yj1(out[, j], lambda[j]) }
    colnames(out) <- paste(colnames(out), round(lambda, 2), sep="^")
    out}  else
    yj1(out, lambda)
  out}

powerTransform <- function(object, ...) UseMethod("powerTransform")

powerTransform.default <- function(object, family="bcPower", ...) {
   y <- object
   if(!inherits(y, "matrix") & !inherits(y, "data.frame")) {
       y <- matrix(y,ncol=1)
       colnames(y) <- c(paste(deparse(substitute(object))))}
   y <- na.omit(y)
   x <- rep(1, dim(y)[1])
   estimateTransform(x, y, NULL, family=family, ...)
   }

powerTransform.lm <- function(object, family="bcPower", ...) {
    mf <- if(is.null(object$model))
            update(object, model=TRUE, method="model.frame")$model
            else object$model
    mt <- attr(mf, "terms")
        y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (is.null(w)) w <- rep(1, dim(mf)[1])
    if (is.empty.model(mt)) {
        x <- matrix(rep(1,dim(mf)[1]), ncol=1) }
    else {
        x <- model.matrix(mt, mf, contrasts)   }
  estimateTransform(x, y, w, family=family, ...)
  }

powerTransform.formula <- function(object, data, subset, weights,
              na.action, family="bcPower", ...) {
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("object", "data", "subset", "weights", "na.action"),
         names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    names(mf)[which(names(mf)=="object")] <- "formula"
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (is.null(w)) w <- rep(1, dim(mf)[1])
    if (is.empty.model(mt)) {
        x <- matrix(rep(1, dim(mf)[1]), ncol=1) }
    else {
        x <- model.matrix(mt, mf)   }
  estimateTransform(x, y, w, family=family, ...)
  }

estimateTransform <- function(X, Y, weights=NULL, family="bcPower", ...) {
  Y <- as.matrix(Y)
  switch(family,
         bcnPower = estimateTransform.bcnPower(X, Y, weights,  ...),
         estimateTransform.default(X, Y, weights, family, ...)
  )
}

# estimateTransform.default is renamed 'estimateTransform
estimateTransform.default <- function(X, Y, weights=NULL,
                    family="bcPower", start=NULL, method="L-BFGS-B", ...) {
  fam <- match.fun(family)
  Y <- as.matrix(Y) # coerces Y to be a matrix.
  X <- as.matrix(X) # coerces X to be a matrix.
  w <- if(is.null(weights)) 1 else sqrt(weights)
  nc <- dim(Y)[2]
  nr <- nrow(Y)
  xqr <- qr(w * X)
  llik <- function(lambda){
    (nr/2)*log(((nr - 1)/nr) *
                 det(var(qr.resid(xqr, w*fam(Y, lambda, j=TRUE, ...)))))
  }
  llik1d <- function(lambda,Y){
    (nr/2)*log(((nr - 1)/nr) * var(qr.resid(xqr, w*fam(Y, lambda, j=TRUE, ...))))
  }
  if (is.null(start)) {
    start <- rep(1, nc)
    for (j in 1:nc){
      res<- suppressWarnings(optimize(
        f = function(lambda) llik1d(lambda,Y[ , j, drop=FALSE]),
        lower=-3, upper=+3))
      start[j] <- res$minimum
    }
  }
  res <- optim(start, llik, hessian=TRUE, method=method,  ...)
  if(res$convergence != 0)
    warning(paste("Convergence failure: return code =", res$convergence))
  res$start<-start
  res$lambda <- res$par
  names(res$lambda) <-
    if (is.null(colnames(Y))) paste("Y", 1:dim(Y)[2], sep="")
  else colnames(Y)
  roundlam <- res$lambda
  stderr <- sqrt(diag(solve(res$hessian)))
  lamL <- roundlam - 1.96 * stderr
  lamU <- roundlam + 1.96 * stderr
  for (val in rev(c(1, 0, -1, .5, .33, -.5, -.33, 2, -2))) {
    sel <- lamL <= val & val <= lamU
    roundlam[sel] <- val
  }
  res$roundlam <- roundlam
  res$invHess <- solve(res$hessian)
  res$llik <- res$value
  res$par <- NULL
  res$family<-family
  res$xqr <- xqr
  res$y <- Y
  res$x <- as.matrix(X)
  res$weights <- weights
  res$family<-family
  class(res) <- "powerTransform"
  res
}

testTransform <- function(object, lambda) UseMethod("testTransform")

testTransform.powerTransform <- function(object, lambda=rep(1, dim(object$y)[2])){
   fam <- match.fun(object$family)
   Y <- cbind(object$y) # coerces Y to be a matrix.
   nc <- dim(Y)[2]
   nr <- nrow(Y)
   lam <- if(length(lambda)==1) rep(lambda, nc) else lambda
   xqr <- object$xqr
   w <- if(is.null(object$weights)) 1 else sqrt(object$weights)
   llik <- function(lambda){
        (nr/2) * log(((nr - 1)/nr) *
             det(var(qr.resid(xqr, w * fam(Y, lam, jacobian.adjusted=TRUE)))))
        }
   LR <- 2 * (llik(lambda) - object$value)
   df <- length(object$lambda)
   pval <- 1-pchisq(LR, df)
   out <- data.frame(LRT=LR, df=df, pval=format.pval(pval))
   rownames(out) <-
     c(paste("LR test, lambda = (",
             paste(round(lam, 2), collapse=" "), ")", sep=""))
   out}

print.powerTransform<-function(x, ...) {
   lambda <- x$lambda
   if (length(lambda) > 1) cat("Estimated transformation parameters \n")
   else cat("Estimated transformation parameter \n")
   print(x$lambda)
   invisible(x)}

summary.powerTransform<-function(object,...){
    one <- 1==length(object$lambda)
    label <- paste(object$family,
       (if(one) "Transformation to Normality" else
                "Transformations to Multinormality"), "\n")
    lambda<-object$lambda
    roundlam <- round(object$roundlam, 2)
    stderr<-sqrt(diag(object$invHess))
    df<-length(lambda)
#    result <- cbind(lambda, roundlam, stderr, lambda - 1.96*stderr, lambda + 1.96*stderr)
    result <- cbind(lambda, roundlam, lambda - 1.96*stderr, lambda + 1.96*stderr)
    rownames(result)<-names(object$lambda)
#    colnames(result)<-c("Est Power", "Rnd Pwr", "Std Err", "Lwr bnd", "Upr Bnd")
    colnames(result)<-c("Est Power", "Rounded Pwr", "Wald Lwr Bnd", "Wald Upr Bnd")
    tests <- testTransform(object, 0)
    tests <- rbind(tests, testTransform(object, 1))
#    if ( !(all(object$roundlam==0) | all(object$roundlam==1) |
#        length(object$roundlam)==1 ))
#           tests <- rbind(tests, testTransform(object, object$roundlam))
    family<-object$family
    out <-  list(label=label, result=result, tests=tests,family=family)
    class(out) <- "summary.powerTransform"
    out
    }

print.summary.powerTransform <- function(x, digits=4, ...) {
    n.trans <- nrow(x$result)
    cat(x$label)
    print(round(x$result, digits))
   if(!is.null(x$family)){
    if(x$family=="bcPower" || x$family=="bcnPower"){
      if (n.trans > 1) cat("\nLikelihood ratio test that transformation parameters are equal to 0\n (all log transformations)\n")
      else cat("\nLikelihood ratio test that transformation parameter is equal to 0\n (log transformation)\n")
      print(x$tests[1,])
      if (n.trans > 1) cat("\nLikelihood ratio test that no transformations are needed\n")
      else cat("\nLikelihood ratio test that no transformation is needed\n")
      print(x$tests[2,])
    }
     if(x$family=="yjPower"){
         if (n.trans > 1) cat("\n Likelihood ratio test that all transformation parameters are equal to 0\n")
         else cat("\n Likelihood ratio test that transformation parameter is equal to 0\n")
       print(x$tests[1,])
     }

   }else{
       if (n.trans > 1) cat("\nLikelihood ratio tests about transformation parameters \n")
       else cat("\nLikelihood ratio test about transformation parameter \n")
     print(x$tests)
   }
}

coef.powerTransform <- function(object, round=FALSE, ...)
  if(round==TRUE) object$roundlam else object$lambda

vcov.powerTransform <- function(object,...) {
  ans <- object$invHess
  rownames(ans) <- names(coef(object))
  colnames(ans) <- names(coef(object))
  ans}









