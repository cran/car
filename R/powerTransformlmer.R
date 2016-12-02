# 2016-07-20:  Added support for power transformations in lmerMod  objects, S. Weisberg
# 2016-08-17:  Allowed fixed gamma or estimated gamma.


# generic functions in powerTransform.R

powerTransform.lmerMod <- function(object, family="bcPower", lambda=c(-3, 3), gamma=NULL, ...) {
  if(family=="skewPower") estimateTransform.skewPowerlmer(object, lambda=lambda, gamma=gamma, ...) else
                          estimateTransform.lmerMod(object, family=family, lambda=lambda, ...)
}

#################################################################################
### estimate transformation methods 
#################################################################################

# lmerMod
estimateTransform.lmerMod <- function(object, family="bcPower", lambda=c(-3, 3), start=NULL, method="L-BFGS-B", ...) {
  data <- model.frame(object)
  y <- (object@resp)$y
  fam <- match.fun(family)
  llik <- function(lambda){
    data$y.lambda <- fam(y, lambda, jacobian.adjusted=TRUE)
    m.lambda <- update(object, y.lambda ~ ., data=data)
    logLik(m.lambda)
  }
  if (is.null(start)) start <- 1
  res<- optimize(f = function(lambda1) llik(lambda1), lower=lambda[1], upper=lambda[2], maximum=TRUE)
# optimize does not give the Hessian, so run optimHess
  res$invHess <- get.invHess(res$maximum, llik, ...)
#  res$hessian <- optimHess(res$maximum, llik, ...)
  res$lambda <- res$maximum
  res$value <- c(res$objective)  
  roundlam <- res$lambda
  stderr <- sqrt(diag(res$invHess))
  lamL <- roundlam - 1.96 * stderr
  lamU <- roundlam + 1.96 * stderr
  for (val in rev(c(1, 0, -1, .5, .33, -.5, -.33, 2, -2))) { 
    sel <- lamL <= val & val <= lamU 
    roundlam[sel] <- val
  }
  res$model <- object
  res$roundlam <- roundlam
  res$family<-family
  class(res) <- c("lmerModpowerTransform", "powerTransform")
  res
}

# skewPowerlmer method
estimateTransform.skewPowerlmer <- function(object, lambda=c(-3, +3), gamma=NULL, ...){
  res <- skewlmermle(object, lambda=lambda, gamma=gamma, ...)
  nc <- 1
  res$lambda <- res$par[1:nc][]
  res$gamma <- res$par[(nc+1):(2*nc)]
  roundlam <- res$lambda
  stderr <- sqrt(diag(res$invHess[1:nc, 1:nc, drop=FALSE]))
  stderr.gam <- sqrt(diag(res$invHess[(nc+1):(2*nc), (nc+1):(2*nc), drop=FALSE]))
  lamL <- roundlam - 1.96 * stderr
  lamU <- roundlam + 1.96 * stderr
  for (val in rev(c(1, 0, -1, .5, .33, -.5, -.33, 2, -2))) { 
    sel <- lamL <= val & val <= lamU 
    roundlam[sel] <- val
  }
  res$roundlam <- roundlam
  res$family <- "skewpowerTransform"
  class(res) <- c("skewpowerTransformlmer", "skewpowerTransform", "powerTransform")
  res  
}

#################################################################################
#  Test Transformation
# in testTransform:  'object' is of class lmerModpowerTransform
#                    'model' will be the lmerMod object
#################################################################################

# lmerMod
testTransform.lmerModpowerTransform <- function(object, lambda=1){
  fam <- match.fun(object$family)
  model <- object$model
  y <- (model@resp)$y
  local.data <- model.frame(model)
  local.data$y.lambda <- fam(y, lambda, jacobian.adjusted=TRUE)
  m.lambda <- update(model, y.lambda ~ ., data=local.data)
  llik <- logLik(m.lambda)
  LR <- c(2 * (object$value - llik))
  df <- 1
  pval <- 1-pchisq(LR, df)
  out <- data.frame(LRT=LR, df=df, pval=pval)
  rownames(out) <- 
    c(paste("LR test, lambda = (",
            paste(round(lambda, 2), collapse=" "), ")", sep=""))
  out}

# skewPower
testTransform.skewpowerTransformlmer <- function(object, lambda=1){
  nc <- 1
  lam <- lambda
  val <- if(is.null(object$fix.gamma))
    skewPowerlmerllikprofile.gamma(object$model, lam)$llik else{
    skewPowerlmerllik(object$model, lam, object$gamma)$llik
    }
  LR <- max(0, -2 * (val - object$llik))
  df <- nc
  pval <- 1-pchisq(LR, df)
  out <- data.frame(LRT=LR, df=df, pval=pval)
  rownames(out) <- 
    c(paste("LR test, lambda = (",
            paste(round(lam, 2), collapse=" "), ")", sep=""))
  out}

###########################################################################
# plot method
###########################################################################
# lmerMod
plot.lmerModpowerTransform <- function(x, z, round=TRUE, plot=pairs, ...){
  cat("plot not supported for mixed models\n") }
# skewPower
plot.skewpowerTransformlmer <- function(x, z, round=TRUE, plot=pairs, ...){
  cat("plot not supported for mixed models\n") }


###########################################################################
##
##  skewPower with lmer objects
##
###########################################################################

## Evaluate skew lmer llik at (lambda, gamma)-----------------------------------
skewPowerlmerllik <- function(object, lambda, gamma) {
  local.data <- model.frame(object)
  Y <- (object@resp)$y
  local.data$y1 <- skewPower(Y, lambda, jacobian.adjusted = TRUE, gamma) 
  mod <- update(object, y1 ~ ., data=local.data)
  list(lambda=lambda, gamma=gamma, llik=logLik(mod))
}

## maximize skewPowerlmerllik for fixed gamma--------------------------------------
skewPowerlmerllikprofile.lambda <- function(object, gamma){
  object <- object
  fn <- function(lam) skewPowerlmerllik(object, lam, gamma)$llik
  f <- optimize(f=fn, interval=c(-3, 3), maximum=TRUE)
  list(lambda=f$maximum, 
       gamma=gamma, 
       llik=f$objective, 
       invHess=solve(-optimHess(f$maximum, fn)))
}

##  maximize skewPowerlmerllik for fixed lambda-------------------------------------
skewPowerlmerllikprofile.gamma <- function(object, lambda){ 
  object <- object
  m <- if(min((object@resp)$y) <= 0) .Machine$double.eps^0.3 else 0L
  fn <- function(gam) skewPowerlmerllik(object, lambda=lambda, gamma=gam)$llik 
  f1 <- optimize(f=fn, interval=c(m, max((object@resp)$y)), maximum=TRUE)
  if (f1$maximum > 10*.Machine$double.eps^0.3){
     hess <- try(optimHess(f1$maximum, fn))
     invHess <- if(class(hess) == "try-error") NA else solve(-hess)} else
     invHess <- NA
  list(lambda=lambda, 
       gamma=f1$maximum, 
       llik=f1$objective, 
       invHess=invHess,
       warn=if(is.na(invHess)) "Estimated gamma is too close to zero; invHess set to NA" else "none")
}

## skewlmerlme  -----------------------------------------------------------
# Use Nelder-Mead simplex to find the mle.  This is slow, but it avoids derivatives that
# cause problems when the mle of gamma is on or near to the boundary of zero.
# method "Nelder-Mead" in optim does not allow box constraints, so I use the 'neldermead'
# function in thneldermead implements the method of Box (1965),
# M. J. Box, "A new method of constrained optimization and a comparison with other methods," 
# Computer J. 8 (1), 42-52 (1965); see also
# http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Nelder-Mead_Simplex
skewlmermle <- function(object, lambda=c(-3, 3), gamma=NULL) {
  if(length(gamma) == 1) skewlmermle1d(object, lambda, gamma) else
    skewlmermle2d(object, lambda, gamma)
}
skewlmermle2d <- function(object, lambda=c(-3, 3), gamma) {
  if(!(requireNamespace("nloptr"))) stop("The 'nloptr' package is needed but is missing")
  # Get a starting value for lambda by removing all negative rows 
  # and doing Box-Cox Power:
  mx <- if(any( (object@resp)$y <= 0)) .Machine$double.eps^0.3 else 0L
  data2 <- model.frame(object)[(object@resp)$y > 0, ]
  if(dim(data2)[2] < 1) stop("Positive responses required to get starting values")
  mod <- update(object, data = data2 )
  lambda.start <- powerTransform(mod, family="bcPower")$lambda
  # Get a starting value for gamma using profiling
  gamma.start <- if(mx == 0L) mx else
    skewPowerlmerllikprofile.gamma(object, lambda.start)$gamma
  pnames <- c("lambda", "gamma")
  fn <- function(param){
    lam <- param[1]
    gam <- param[2]
    skewPowerlmerllik(object, lam, gam)$llik
  } 
  res <- nloptr::neldermead(c(lambda.start, gamma.start), function(param) -fn(param),  
                            lower=c(lambda[1], mx), 
                            upper=c(lambda[2], +Inf)) 
  fit <- NULL
  # check for convergence
  fit$message <- if(res$convergence==0)  
    "Normal convergence" else 
      paste("Convergence likely at a boundary point, convergence code:", res$convergence)
  # change sign of res$value
  fit$par <- res$par
  names(fit$par) <- pnames
  fit$llik <- -res$value
  # Compute the inverse Hessian if possible using optimHess
  hess <- try(optimHess(res$par, fn), silent=TRUE)
  if(class(hess) == "try-error"){ 
    fit$warn <- "estimate of gamma is very close to 0 -- Hessian computed for lambda only"
    fit$invHess <- matrix(NA, nrow=2, ncol=2)
    fit$invHess[1, 1] <- skewPowerlmerllikprofile.lambda(object, fit$par[2])$invHess} 
  else{
    fit$invHess <- solve(-hess)
    fit$fix.gamma <- NULL
    fit$warn <- "none"
  }
  # add names
  dimnames(fit$invHess) <- list(pnames, pnames)
  fit$ylabs <- "Y" 
  fit$model <- object
  fit
}

# Skew power with gamma fixed, not estimated
skewlmermle1d <- function(object, lambda, gamma) {
  pnames <- c("lambda", "gamma")
  res <- skewPowerlmerllikprofile.lambda(object, gamma)
  fit <- NULL
  fit$par <- c(lambda=res$lambda, gamma=gamma)
  fit$llik <- res$llik
  fit$invHess <- matrix(c(res$invHess, NA, NA, NA), nrow=2, ncol=2, dimnames = list(pnames, pnames))
  names(fit$par) <- pnames
  fit$ylabs <- "Y" 
  fit$model <- object
  fit$fix.gamma <- gamma
  fit
}


# methods for skewpowerTransform objects

summary.skewpowerTransformlmer<-function(object,...){
  nc <- length(object$lambda)
  label <- "Skew Power transformation to Normality, lmer fit\n\n" 
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
  out <-  list(label=label, result=result, result.gamma=result.gamma, tests=tests)
  class(out) <- "summary.skewpowerTransform"
  out
}


## ------------------------------------------------------------------------
contour.skewpowerTransformlmer <- function(x, ksds=4, levels=c(.5, .95, .99, .999), 
                                           main="Skew Power Log-likelihood", ...){
  object <- x
  if(dim(object$y)[2] != 1L) stop("This function is for univariate Y only")
  q <- object$value - qchisq(levels, 2)/2
  se <- sqrt(diag(object$invHess))
  center <- c(object$lambda, object$gamma)
  x1 <- seq(object$lambda - ksds*se[1], object$lambda + ksds*se[1], length=100) 
  y <- seq(max(.01, object$gamma - ksds*se[2]), object$gamma + ksds*se[2], length=100)
  z <- matrix(0, nrow=length(x1), ncol=length(y))
  for (i in 1:length(x1)){
    for (j in 1:length(y)){
      z[i,j] <- skewPowerlmerllik(object, x1[i], y[j])$llik
    }
  }
  contour(x1, y, z, xlab=expression(lambda), ylab=expression(gamma), main=main, 
          nlevels=length(levels), levels=q, ...)
  points(center[1], center[2], pch=16, cex=1.25)
  text(center[1], center[2], as.character(round(object$value, 2)), pos=4, cex=.75)
}

###########################################################################
### boxCox method
###########################################################################
# Note:
# boxCox.lm works correctly with the skewPower family, so there is no separate method
# lmerMod method also works with skewPower
boxCox.lmerMod <- function(object,
                           lambda = seq(-2, 2, 1/10), plotit = TRUE, 
                           interp = plotit, eps = 1/50,
                           xlab=NULL, ylab=NULL, 
                           family="bcPower", 
                           param=c("lambda", "gamma"), gamma=NULL, grid=TRUE, ...)
{ 
  param <- match.arg(param)
  ylab <- if(is.null(ylab)){if(family != "skewPower") "log-likelihood" else{
    if(param=="gamma") {expression(max(logL[gamma](lambda,gamma)))} else 
    {expression(max[lambda](logL(lambda, gamma)))}}} else ylab
  xlab <- if(is.null(xlab)){if(param == "lambda") expression(lambda) else expression(gamma)} else xlab
  fam <- match.fun(family)
  fdata <- model.frame(object)
  y <- (object@resp)$y
  fam <- match.fun(family)
  xl <- loglik <- if(family != "skewPower") as.vector(lambda) else {
    if(param == "lambda") as.vector(lambda) else {
      # if argument gamma is non-null, use it for the range for gamma.  
      # if gamma is null then use the range of the mle plus or minus 3 ses
      if(!is.null(gamma)) as.vector(gamma) else{ 
        p1 <- powerTransform(object, family="skewPower")
        gam <- p1$gamma
        se <- sqrt(vcov(p1)[2,2])
        seq(max(.01, gam - 3*se), gam + 3*se, length=100)
      }
    }
  } 
  m <- length(xl)
  loglik <- rep(0, m)
  if(family != "skewPower"){
    for (i in 1L:m) {
      fdata$y.lambda <- fam(y, xl[i], jacobian.adjusted=TRUE)
      m.lambda <- update(object, y.lambda ~ ., data=fdata)
      loglik[i] <- c(logLik(m.lambda))
    }}  else{ 
      for (i in 1L:m) { 
        loglik[i] <- if(param == "gamma")
            skewPowerlmerllikprofile.lambda(object, xl[i])$llik else
            skewPowerlmerllikprofile.gamma(object, xl[i])$llik 
      }
    }
  if (interp) {
    sp <- spline(xl, loglik, n = 100)
    xl <- sp$x
    loglik <- sp$y
    m <- length(xl)
  }
  if (plotit) {
    mx <- (1L:m)[loglik == max(loglik)][1L]
    Lmax <- loglik[mx]
    lim <- Lmax - qchisq(19/20, 1)/2
    plot(xl, loglik, xlab = xlab, ylab = ylab, type = "n",
         ylim = range(loglik, lim), ...)
    if(grid){
      grid(lty=1, equilogs=FALSE)
      box()}
    lines(xl, loglik)
    plims <- par("usr")
    abline(h = lim, lty = 2)
    y0 <- plims[3L]
    scal <- (1/10 * (plims[4L] - y0))/par("pin")[2L]
    scx <- (1/10 * (plims[2L] - plims[1L]))/par("pin")[1L]
    text(xl[1L] + scx, lim + scal, " 95%")
    la <- xl[mx]
    if (mx > 1 && mx < m)
      segments(la, y0, la, Lmax, lty = 2)
    ind <- range((1L:m)[loglik > lim])
    if (loglik[1L] < lim) {
      i <- ind[1L]
      x <- xl[i - 1] + ((lim - loglik[i - 1]) * (xl[i] -
                                                   xl[i - 1]))/(loglik[i] - loglik[i - 1])
      segments(x, y0, x, lim, lty = 2)
    }
    if (loglik[m] < lim) {
      i <- ind[2L] + 1
      x <- xl[i - 1] + ((lim - loglik[i - 1]) * (xl[i] -
                                                   xl[i - 1]))/(loglik[i] - loglik[i - 1])
      segments(x, y0, x, lim, lty = 2)
    }
  }
  invisible(list(x = xl, y = loglik))
}



