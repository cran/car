# 2015-08-26:  Modified by S. Weisberg to add support for bcn power transformations.
# 2017-05-11:  Added boxCox2d, renamed verssion of contour.powerTransform
# 2017-05-11:  Bug fixes in boxCox.formula with arugment passing to other methods

boxCox <- function(object,...) UseMethod("boxCox")

#  New arguments:  param, and gamma
boxCox.formula <- function (object, lambda = seq(-2, 2, 1/10), plotit=TRUE,
          family="bcPower", param=c("lambda", "gamma"), gamma=NULL, grid=TRUE,...)
{
  m <- length(lambda)
  object <- lm(object, y = TRUE, qr = TRUE)
  result <- NextMethod()
  if (plotit)
    invisible(result)
  else result
}

boxCox.lm <- function (object, lambda = seq(-2, 2, 1/10),
                       plotit = TRUE, ...)
{
  m <- length(lambda)
  if (is.null(object$y) || is.null(object$qr))
    object <- update(object, y = TRUE, qr = TRUE)
  result <- NextMethod()
  if (plotit)
    invisible(result)
  else result
}

boxCox.default <- function(object,
                           lambda = seq(-2, 2, 1/10), plotit = TRUE,
                           interp = plotit, eps = 1/50,
                           xlab=NULL, ylab=NULL,
                           family="bcPower",
                           param=c("lambda", "gamma"), gamma=NULL, grid=TRUE, ...)
{
  if(class(object)[1] == "mlm")  stop("This function is for univariate response only")
  param <- match.arg(param)
  ylab <- if(is.null(ylab)){if(family != "bcnPower") "log-likelihood" else{
    if(param=="gamma") {expression(max(logL[gamma](lambda,gamma)))} else
    {expression(max[lambda](logL(lambda, gamma)))}}} else ylab
  xlab <- if(is.null(xlab)){if(param == "lambda") expression(lambda) else expression(gamma)} else xlab
  fam <- match.fun(family)
  if (is.null(object$y) || is.null(object$qr))
    stop(paste(deparse(substitute(object)),
               "does not have both 'qr' and 'y' components"))
  y <- object$y
  n <- length(y)
  xqr <- object$qr
  xl <- loglik <- if(family != "bcnPower") as.vector(lambda) else {
    if(param == "lambda") as.vector(lambda) else {
# if argument gamma is non-null, use it for the range for gamma.
# if gamma is null then use the range of the mle plus or minus 3 ses
      if(!is.null(gamma)) as.vector(gamma) else{
        p1 <- powerTransform(object, family="bcnPower")
        gam <- p1$gamma
        se <- sd(y)  # arbitrary scaling factor
        seq(max(.01, gam - 3*se), gam + 3*se, length=100)
      }
    }
  }
  m <- length(xl)
  if(family != "bcnPower"){
    for (i in 1L:m) {
      yt <- fam(y,xl[i],j=TRUE)
      loglik[i] <- -n/2 * log(sum(qr.resid(xqr, yt)^2))
    }}  else{
      lambda.1d <- function(gamma){
        fn <- function(lam) bcnPowerllik(NULL, y, NULL, lambda=lam, gamma=gamma, xqr=xqr)$llik
        f <- optimize(f=fn, interval=c(-3, 3), maximum=TRUE)
        f$objective
      }
      gamma.1d <- function(lambda){
        fn <- function(gam) bcnPowerllik(NULL, y, NULL, lambda=lambda, gamma=gam, xqr=xqr)$llik
        f <- optimize(f=fn, interval=c(0.01, max(y)), maximum=TRUE)
        f$objective
      }
      for (i in 1L:m) {
        loglik[i] <- if(param == "lambda")
          gamma.1d(loglik[i]) else
           lambda.1d(loglik[i])
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
         ylim = range(loglik, lim))
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
  list(x = xl, y = loglik)
}


##########
boxCox2d <- function(x, ksds=4, levels=c(.5, .95, .99, .999),
                                      main="bcnPower Log-likelihood", grid=TRUE, ...){
  if(class(x)[1] != "bcnPowerTransform")
    stop("Error--first argument must be a bcnPower transformation")
  object <- x
  if(dim(object$y)[2] != 1L) stop("This function is for univariate Y only")
  q <- object$llik - qchisq(levels, 2)/2
  se <- sqrt(diag(object$invHess))
  center <- c(object$lambda, object$gamma)
  x1 <- seq(object$lambda - ksds*se[1], object$lambda + ksds*se[1], length=100)
  y <- seq(max(.01, object$gamma - ksds*se[2]), object$gamma + ksds*se[2], length=100)
  z <- matrix(0, nrow=length(x1), ncol=length(y))
  for (i in 1:length(x1)){
    for (j in 1:length(y)){
      z[i,j] <- bcnPowerllik(NULL, object$y, object$weights, x1[i], y[j], xqr=object$xqr)$llik
    }
  }
  contour(x1, y, z, xlab=expression(lambda), ylab=expression(gamma), main=main,
          nlevels=length(levels), levels=q, ...)

  points(center[1], center[2], pch=16, cex=1.25)
  text(center[1], center[2], as.character(round(object$llik, 2)), pos=4, cex=.75)
  if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
}
