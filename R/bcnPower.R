# 05-02-2017:  bcnPower family, replacing skewPower.  S. Weisberg
# 2017-05-18: Changed summary.powerTransform; deleted invalid test; added roundlam to output
# 2017-12-19: Deleted plot method
# 2017-12-19: Improved handling of gamma small case, still not great for the
#             multivariate extenstion.  Works for lm and lmer
# 2017-12-25: bug fix with multivariace bcnPower
# 2019-03-07: bug fix in estimateTransform.bcnPowerlmer, thanks to wouter@zoology.ubc.ca
# 2019-11-14,15: change class(x) == "y" to inherits(x, "y") and likewise for !=

bcnPower <- function(U, lambda, jacobian.adjusted=FALSE, gamma) {
  if(is.matrix(U)){
    if(dim(U)[2] != length(lambda) | dim(U)[2] != length(gamma))
      stop("gamma and lambda must have length equal to number of columns in U")
  } else {
    if(length(gamma) != 1 | length(lambda) != 1)
      stop("gamma and lambda must be length 1")
  }
  if(any(gamma < 0)) stop("gamma must be >= 0")
  hc1 <- function(U, lambda, gamma){
    if(abs(gamma) <= 1.e-10 & any(U[!is.na(U)] <= 0))
      stop("First argument must be strictly positive if gamma = 0.")
    s <- sqrt(U^2 + gamma^2)
    z <- if (abs(lambda) <= 1.e-10)
      log(.5*(U + s)) else ((.5*(U + s))^lambda - 1)/lambda
    if (jacobian.adjusted == TRUE) {
      Jn <- (.5^lambda) *
        (exp((lambda - 1) * mean(log(U + s), na.rm=TRUE))) *
        (exp(mean(log(1 + U/s), na.rm=TRUE)))
      z <- z/Jn}
    z
  }
  out <- U
  out <- if(is.matrix(out) | is.data.frame(out)){
    if(is.null(colnames(out)))
      colnames(out) <- paste("Z", 1:dim(out)[2], sep="")
    for (j in 1:ncol(out)) {out[, j] <- hc1(out[, j], lambda[j], gamma[j]) }
    colnames(out) <- paste(colnames(out), "(",
                           round(lambda, 2), ",",round(gamma, 1),")", sep="")
    #    colnames(out) <- paste(colnames(out), round(lambda, 2), sep="^")
    out}  else
      hc1(out, lambda, gamma)
  out}

bcnPowerInverse <- function(z, lambda, gamma){
  q <- if(abs(lambda) < 1.e-7) 2 * exp(z) else 2 * (lambda*z + 1)^(1/lambda)
  (q^2 - gamma^2)/(2 * q)
}

###############################################################################
# estimateTransform and methods
#

# multivariate box-cox with negatives starting values, given X and Y
bcn.sv <- function(X, Y, weights, itmax=100, conv=.0001, verbose=FALSE,
                   start=TRUE, gamma.min=.1){
  Y <- as.matrix(Y)
  d <- dim(Y)[2]
  if(d > 1) stop("bcn.sv requires a univariate response")
  lambda.1d <- function(Y, weights, lambda, gamma, xqr){
    fn <- function(lam) bcnPowerllik(NULL, Y, weights, lambda=lam,
                                     gamma=gamma, xqr=xqr)$llik
    f <- optimize(f=fn, interval=c(-3, 3), maximum=TRUE)
    list(lambda=f$maximum, gamma=gamma, llik=f$objective)
  }
  gamma.1d <- function(Y, weights, lambda, gamma, xqr){
    fn1 <- function(gam) bcnPowerllik(NULL, Y, weights, lambda=lambda,
                                     gamma=gam, xqr=xqr)$llik
    f <- optimize(f=fn1, interval=c(0.01, max(Y)), maximum=TRUE)
    list(lambda=lambda, gamma=f$maximum, llik=f$objective)
  }
  # get qr decomposition
  w <- if(is.null(weights)) 1 else sqrt(weights)
  xqr <- qr(w * as.matrix(X))
  # get starting value for gamma
  gamma <- if(min(Y) <= 0) max(min(Y[Y>0]), 5*gamma.min) else 0
  res <- lambda.1d(Y, weights, lambda=1, gamma=gamma, xqr)
  res <- gamma.1d(Y, weights, lambda=res$lambda, gamma=res$gamma, xqr)
  # set iteration counter
  i <- 0
  crit <- 1
  gamma.ok <- TRUE
  while( (crit > conv) & (i < itmax) & gamma.ok) {
    i <- i+1
    last.value <- res
    res <- lambda.1d(Y, weights, res$lambda, res$gamma, xqr)
    res <- gamma.1d(Y, weights, res$lambda, res$gamma, xqr)
    if(res$gamma < 1.5 * gamma.min){
      gamma.ok <- FALSE
      res <- lambda.1d(Y, weights, res$lambda, gamma.min, xqr)
    }
    crit <- (res$llik - last.value$llik)/abs(res$llik)
    if(verbose)
      print(data.frame(Iter=i, gamma=res$gamma,
            lambda=res$gamma, llik=res$llik, crit=crit))
  }
  if(i==itmax & conv > crit)
    warning(paste("No convergence in", itmax, "iterations, criterion =", crit, collapse=" "))
#  if(!gamma.ok) warning(paste("gamma too close to zero, set to", gamma.min, collapse=" "))
  if(start == TRUE)  return(c(res, gamma.estimated=gamma.ok)) else {
# compute the Hessian -- depends on gamma.ok
  if(gamma.ok){
    fn2 <- function(param){
      lam <- param[1]
      gam <- param[2]
      bcnPowerllik(NULL, Y, weights, lam, gam, xqr=xqr)$llik
    }
    hess <- optimHess(c(res$lambda, res$gamma), fn2)
    res$invHess <- solve(-hess)} else{
# gamma.ok == FALSE
      fn3 <- function(lam){
        lam
        bcnPowerllik(NULL, Y, weights, lam, gamma.min, xqr=xqr)$llik
      }
      hess <- optimHess(res$lambda, fn3)
      res$invHess <- matrix(c(-1/hess, NA, NA, NA), ncol=2)}
# end computing of invHess
    rownames(res$invHess) <- colnames(res$invHess) <- c("lambda", "gamma")
    roundlam <- res$lambda
    stderr <- sqrt(diag(res$invHess[1, 1, drop=FALSE]))
    stderr.gam <- sqrt(diag(res$invHess[2, 2, drop=FALSE]))
    lamL <- roundlam - 1.96 * stderr
    lamU <- roundlam + 1.96 * stderr
    for (val in rev(c(1, 0, -1, .5, .33, -.5, -.33, 2, -2))) {
      sel <- lamL <= val & val <= lamU
      roundlam[sel] <- val
    }
    res$roundlam <- roundlam
    res$ylabs <-
      if (is.null(colnames(Y))) paste("Y", 1:dim(as.matrix(Y))[2], sep="") else colnames(Y)
    res$xqr <- xqr
    res$y <- as.matrix(Y)
    res$x <- as.matrix(X)
    res$weights <- weights
    res$family <- "bcnPowerTransform"
    res$y
    res$gamma.estimated <- gamma.ok
    class(res) <- c("bcnPowerTransform", "powerTransform")
    res}
}

estimateTransform.bcnPower <- function(X, Y, weights,
          itmax=100, conv=.0001, verbose=FALSE, gamma.min=.1){
  d <- dim(as.matrix(Y))[2]
  skf.lambda <- function(Y, weights, lambda, gamma, xqr){
    fn3a <- function(lam) bcnPowerllik(NULL, Y, weights, lambda=lam,
                                     gamma=gamma, xqr=xqr)$llik
    f <- optim(par=lambda, fn=fn3a, method="L-BFGS-B",
                 lower=rep(-3, d), upper=rep(3, d),
                 control=list(fnscale=-1))
    list(lambda=f$par, gamma=gamma, llik=f$value, conv=f$convergence, message=f$message)
  }
  skf.gamma <- function(Y, weights, lambda, gamma, xqr){
    fn3b <- function(gam) bcnPowerllik(NULL, Y, weights, lambda=lambda,
                                     gamma=gam, xqr=xqr)$llik
    f <- optim(par=gamma, fn=fn3b, method="L-BFGS-B",
                 lower=rep(gamma.min, d),
                 upper=rep(Inf, d),
                 control=list(fnscale=-1))
      list(lambda=lambda, gamma=f$par, llik=f$value,
           conv=f$convergence, message=f$message)
  }
# get qr decomposition once
  w <- if(is.null(weights)) 1 else sqrt(weights)
  xqr <- qr(w * as.matrix(X))
# if d = 1 call bcn.sv and return, else call bcn.sv to get starting values.
  if(d == 1) bcn.sv(X, Y, weights, start=FALSE) else{
# The rest of this code is for the multivariate case
# get starting values for gamma
  sv <- apply(Y, 2, function(y) unlist(bcn.sv(X, y, weights, start=TRUE)))
  res <- as.list(as.data.frame(t(sv))) # output to a list
# gamma.estimated converted to numeric, so fixup
  res$gamma.estimated <- ifelse(res$gamma.estimated==1, TRUE, FALSE)
  res$llik <- -Inf
# set iteration counter
  i <- 0
  crit <- 1
# iterate
  while( (crit > conv) & (i < itmax)) {
    i <- i+1
    last.value <- res
    res <- skf.gamma (Y, weights, res$lambda, res$gamma, xqr)
    res <- skf.lambda(Y, weights, res$lambda, res$gamma, xqr)
    crit <- (res$llik - last.value$llik)/abs(res$llik)
    if(verbose)
      print(paste("Iter:", i, "llik=", res$llik, "Crit:", crit, collapse=" "))
  }
  if(itmax == 1) warning("One iteration only, results assume responses are uncorrelated")
#  if(i==itmax & conv > crit)
#    warning(paste("No convergence in", itmax, "iterations, criterion =", crit, collapse=" "))
  fn4 <- function(param){
    lam <- param[1:d]
    gam <- param[(d+1):(2*d)]
    bcnPowerllik(NULL, Y, weights, lam, gam, xqr=xqr)$llik
  }
# check gamma
  gamma.ok <- ifelse(res$gamma > 1.5*gamma.min, TRUE, FALSE)
  res$gamma[!gamma.ok] <- gamma.min
  if(all(gamma.ok)){
     hess <- try(optimHess(c(res$lambda, res$gamma), fn4))
     res$invHess <- if(inherits(hess, "try-error")) NA else solve(-hess)
  } else {
    fn4a <- function(lam) fn4(c(lam, res$gamma))
    hess <- try(optimHess(res$lambda, fn4a)) # hessian for lambda only
    res$invHess <- matrix(NA, nrow=2*d, ncol=2*d)
    res$invHess[1:d, 1:d] <- solve(-hess)
    }
  roundlam <- res$lambda
  stderr <- sqrt(diag(res$invHess[1:d, 1:d, drop=FALSE]))
  stderr.gam <- sqrt(diag(res$invHess[(d+1):(2*d), (d+1):(2*d), drop=FALSE]))
  lamL <- roundlam - 1.96 * stderr
  lamU <- roundlam + 1.96 * stderr
  for (val in rev(c(1, 0, -1, .5, .33, -.5, -.33, 2, -2))) {
    sel <- lamL <= val & val <= lamU
    roundlam[sel] <- val
  }
  res$roundlam <- roundlam
  res$ylabs <-
    if (is.null(colnames(Y))) paste("Y", 1:d, sep="") else colnames(Y)
  invHesslabels <- c(paste(res$ylabs, "lambda", sep=":"),
                     paste(res$ylabs, "gamma", sep=":"))
  if (!inherits(hess, "try-error"))
        rownames(res$invHess) <- colnames(res$invHess) <- invHesslabels
  res$xqr <- xqr
  res$y <- as.matrix(Y)
  res$x <- as.matrix(X)
  res$weights <- weights
  res$family <- "bcnPowerTransform"
  res$y
  class(res) <- c("bcnPowerTransform", "powerTransform")
  res$gamma.estimated <- gamma.ok
  res
}}

#############################################################################
## The log-likelihood function assuming a normal target
## Evaluate bcnPower llik at (lambda, gamma)-----------------------------------
bcnPowerllik <- function(X, Y, weights=NULL, lambda, gamma, xqr=NULL) {
  Y <- as.matrix(Y) # coerces Y to be a matrix.
  w <- if(is.null(weights)) 1 else sqrt(weights)
  xqr <- if(is.null(xqr)){qr(w * as.matrix(X))} else xqr
  nr <- nrow(Y)
  f <- -(nr/2)*log(((nr - 1)/nr) *
      det(as.matrix(var(qr.resid(xqr, w * bcnPower(Y, lambda,
                          jacobian.adjusted=TRUE, gamma=gamma))))))
  list(lambda=lambda, gamma=gamma, llik=f)
}

###############################################################################
# testTransform
testTransform.bcnPowerTransform <- function(object, lambda=rep(1, dim(object$y)[2])){
  d <- length(object$lambda)
  lam <- if(length(lambda)==1) rep(lambda, d) else lambda
  skf.gamma <- function(Y, weights, lambda, gamma, xqr){
    fn5 <- function(gam) bcnPowerllik(NULL, Y, weights, lambda=lam,
                                     gamma=gamma, xqr=xqr)$llik
    f <- optim(par=gamma, fn=fn5, method="L-BFGS-B",
               lower=rep(.Machine$double.eps^0.25, d),
               upper=rep(Inf, d),
               control=list(fnscale=-1))
    list(lambda=lambda, gamma=f$par, llik=f$value, conv=f$convergence,
         message=f$message)
  }
  val <- skf.gamma(object$y, object$weights, lam,
         gamma=object$gamma, xqr=object$xqr)$llik
  LR <- max(0, -2 * (val - object$llik))
  df <- d
  pval <- 1-pchisq(LR, df)
  out <- data.frame(LRT=LR, df=df, pval=pval)
  rownames(out) <-
    c(paste("LR test, lambda = (",
            paste(round(lam, 2), collapse=" "), ")", sep=""))
  out}

print.bcnPowerTransform<-function(x, ...) {
  cat("Estimated transformation power, lambda\n")
  print(x$lambda)
# temporary code
  if(is.null(x$gamma.estimated)) x$gamma.estimated=TRUE
  if(any(x$gamma.estimated)){
    cat("\nEstimated location, gamma\n")} else{
    cat("\nLocation gamma was fixed at its lower bound\n")}
  print(x$gamma)
  invisible(x)}

summary.bcnPowerTransform <- function(object, ...){
  nc <- length(object$lambda)
  label <- paste(if(nc==1) "bcnPower transformation to Normality" else
                   "bcnPower transformation to Multinormality", "\n")
  lambda <- object$lambda
  roundlam <- round(object$roundlam, 3)
  gamma <- object$gamma
  stderr <- sqrt(diag(object$invHess))
  stderr.gamma <- stderr[(nc+1):(2*nc)]
  stderr <- stderr[1:nc]
  result <- cbind(lambda, roundlam, lambda - 1.96*stderr, lambda + 1.96*stderr)
  result.gamma <- cbind(gamma, stderr.gamma, pmax(gamma - 1.96*stderr.gamma, 0), gamma + 1.96*stderr.gamma)
  rownames(result) <- rownames(result.gamma) <- object$ylabs
  colnames(result) <- c("Est Power", "Rounded Pwr", "Wald Lwr Bnd", "Wald Upr Bnd")
  colnames(result.gamma) <-
    c("Est gamma", "Std Err.", "Wald Lower Bound", "Wald Upper Bound")
  tests <- testTransform(object, 0)
  tests <- rbind(tests, testTransform(object, 1))
#  if ( !(all(object$roundlam==0) | all(object$roundlam==1) |
#           length(object$roundlam)==1 | all(object$roundlam == object$lambda)))
#    tests <- rbind(tests, testTransform(object, object$roundlam))
  out <-  list(label=label, result=result, result.gamma=result.gamma,
               tests=tests, gamma.estimated=object$gamma.estimated)
  if(is.null(out$gamma.estimated)) out$gamma.estimated <- TRUE
  class(out) <- "summary.bcnPowerTransform"
  out
}

print.summary.bcnPowerTransform <- function(x,digits=4, ...) {
  cat(x$label)
  cat("\nEstimated power, lambda\n")
  print(round(x$result, digits))
  if(any(x$gamma.estimated)){
     cat("\nEstimated location, gamma\n")} else{
     cat("\nLocation gamma was fixed at its lower bound\n")}
  print(round(x$result.gamma, digits))
  cat("\nLikelihood ratio tests about transformation parameters\n")
  print(x$tests)
  if(any(x$result.gamma[,1] < 1.e-5)) warning(
    "When gamma is zero, transformation family is the Box-Cox Power family")
}

coef.bcnPowerTransform <- function(object, param=c("both", "lambda", "gamma"), round=FALSE, ...){
  param <- match.arg(param)
  co <- cbind(if(round==TRUE) object$roundlam else object$lambda, object$gamma)
  dimnames(co) <- list(object$ylabs, c("lambda", "gamma"))
  switch(param, lambda = co[, 1], gamma=co[, 2], both= co)
}

vcov.bcnPowerTransform <- function(object, param=c("both", "lambda", "gamma"), ...) {
  param <- match.arg(param)
  nc <- length(object$lambda)
  switch(param, lambda=object$invHess[1:nc, 1:nc], gamma=object$invHess[(nc+1):(2*nc), (nc+1):(2*nc)],
         both=object$invHess)
}





##########################################################################################
#  bcnPower for lmer models
#  Modified 12/19/2017 to handle gamma-at-the boundary gracefully
estimateTransform.bcnPowerlmer <- function(object, verbose=FALSE,
            conv=.001, itmax=100, gamma.min=.1, ...) {
  data <- model.frame(object)
  y <- (object@resp)$y
  lambda.1d <- function(lambda, gamma){
    fn6 <- function(lam){
      data$y1 <- bcnPower(y, lambda=lam, jacobian.adjusted=TRUE, gamma)
      logLik(update(object, y1 ~ ., data=data))}
    f <- optimize(f=fn6, interval=c(-3, 3), maximum=TRUE)
    list(lambda=f$maximum, gamma=gamma, llik=f$objective)
  }
  gamma.1d <- function(lambda=lambda, gamma=gamma){
    fn7 <- function(gam){
      data$y1 <- bcnPower(y, lambda, jacobian.adjusted=TRUE, gamma=gam)
      logLik(update(object, y1 ~ ., data=data))}
    f <- optimize(f=fn7, interval=c(.5*gamma.min, max(y)), maximum=TRUE)
    list(lambda=lambda, gamma=f$maximum, llik=f$objective)
  }
  # starting values for lambda, gamma
  lambda <- gamma <- 1
  gamma.ok <- TRUE
  res <- lambda.1d(lambda, gamma)
  res <- gamma.1d(res$lambda, res$gamma)
  if(res$gamma < 1.5 * gamma.min){
    gamma.ok <- FALSE
    res <- lambda.1d(res$lambda, gamma.min)
  } else{
  # iteration is needed only if gamma is not on the boundary
  # set iteration counter
    i <- 0
    crit <- 1
    while( (crit > conv) & (i < itmax) & gamma.ok) {
      i <- i+1
      last.value <- res
      res <- lambda.1d(res$lambda, res$gamma)
      res <- gamma.1d(res$lambda, res$gamma)
      if(res$gamma < 1.5 * gamma.min){
        gamma.ok <- FALSE
        res <- lambda.1d(res$lambda, gamma.min)
      }
      crit <- (res$llik - last.value$llik)/abs(res$llik)
      if(verbose)
      print(data.frame(Iter=i, gamma=res$gamma, lambda=res$lambda,
                       llik=res$llik, crit=crit))
    }
    if(i==itmax & conv > crit)
      warning(paste("No convergence in", itmax, "iterations, criterion =",
                  crit, collapse=" "))
  }
#  if(!gamma.ok) warning(paste("gamma too close to zero, set to",gamma.min, collapse=" "))
  # optimize does not give the Hessian, so run optimHess
  if(gamma.ok){
    llikfn <- function(par){
       data$y1 <- bcnPower(y, par[1], jacobian.adjusted=TRUE, par[2])
       mf <- update(object, y1 ~ ., data=data)
       logLik(mf)
    }
    res$invHess <- solve(-optimHess(unlist(res[1:2]), llikfn))
    if(any(diag(res$invHess) < 0)) res$invHess <- matrix(NA, nrow=2, ncol=2)
    } else
    {
    llikfn1 <- function(lam){
       data$y1 <- bcnPower(y, lambda=lam, jacobian.adjusted=TRUE, gamma=res$gamma)
       logLik(update(object, y1 ~ ., data=data))}
    v1 <- -1/optimHess(res$lambda, llikfn1)
    res$invHess <- matrix(c(v1, NA, NA, NA), ncol=2)
    }
  roundlam <- res$lambda
  stderr <- sqrt(res$invHess[1,1])
  lamL <- roundlam - 1.96 * stderr
  lamU <- roundlam + 1.96 * stderr
  for (val in rev(c(1, 0, -1, .5, .33, -.5, -.33, 2, -2))) {
    sel <- lamL <= val & val <= lamU
    roundlam[sel] <- val
  }
  res$model <- object
  res$roundlam <- roundlam
  res$family<-family
  res$gamma.estimated <- gamma.ok
  class(res) <- c("bcnPowerTransformlmer", "bcnPowerTransform")
  res
}

testTransform.bcnPowerTransformlmer <- function(object, lambda=1){
  nc <- 1
  lam <- lambda
  mod <- object$model
  data <- model.frame(mod)
  data$.y <- mod@resp$y
  gamma.1d <- function(mod, lambda=lambda, gamma=gamma){
    fn <- function(gam){
      data$.y1 <- bcnPower(data$.y, lambda, jacobian.adjusted=TRUE, gamma=gam)
      logLik(update(mod, .y1 ~ ., data=data))}
    f <- optimize(f=fn, interval=c(1.e-5, max(data$.y)), maximum=TRUE)
    list(lambda=lambda, gamma=f$maximum, llik=f$objective)
  }
  val <- gamma.1d(object$model, lambda, object$gamma)$llik
  LR <- max(0, -2 * (val - object$llik))
  df <- nc
  pval <- 1-pchisq(LR, df)
  out <- data.frame(LRT=LR, df=df, pval=pval)
  rownames(out) <-
    c(paste("LR test, lambda = (",
            paste(round(lam, 2), collapse=" "), ")", sep=""))
  out}

summary.bcnPowerTransformlmer<-function(object,...){
  nc <- length(object$lambda)
  label <- "bcn - Box-Cox Power transformation to Normality\nallowing for negative values, lmer fit\n"
  lambda <- object$lambda
  gamma <- object$gamma
  stderr <- sqrt(diag(object$invHess))
  stderr.gamma <- stderr[(nc+1):(2*nc)]
  stderr <- stderr[1:nc]
  result <- cbind(lambda, stderr, lambda - 1.96*stderr, lambda + 1.96*stderr)
  result.gamma <- cbind(gamma, stderr.gamma, pmax(gamma - 1.96*stderr.gamma, 0), gamma + 1.96*stderr.gamma)
  rownames(result) <- rownames(result.gamma) <- object$ylabs
  colnames(result) <- colnames(result.gamma) <-
    c("Est.Power", "Std.Err.", "Wald Lower Bound", "Wald Upper Bound")
  colnames(result.gamma) <-
    c("Est.gamma", "Std.Err.", "Wald Lower Bound", "Wald Upper Bound")
  tests <- testTransform(object, 0)
  tests <- rbind(tests, testTransform(object, 1))
  if ( !(all(object$roundlam==0) | all(object$roundlam==1) |
         length(object$roundlam)==1 ))
    tests <- rbind(tests, testTransform(object, object$roundlam))
  out <-  list(label=label, result=result, result.gamma=result.gamma,
               gamma.estimated=object$gamma.estimated,tests=tests)
  class(out) <- "summary.bcnPowerTransform"
  out
}




