# 2016-07-20:  Added support for power transformations in lmerMod  objects, S. Weisberg
# 2016-05-02:  Moved (working) cosde for bncPower family to bcnPower.R
# 2017-12-19:  Modified estimateTransform to handle gamma \approx 0 gracefully.
# 2017-12-19:  added error for 'I' terms in formulas

# generic functions in powerTransform.R

powerTransform.lmerMod <- function(object, family="bcPower", ...) {
  if(family=="bcnPower") estimateTransform.bcnPowerlmer(object,  ...) else
    estimateTransform.lmerMod(object, family=family, ...)
}

#################################################################################
### estimate transformation methods
#################################################################################

# lmerMod
estimateTransform.lmerMod <- function(object, family="bcPower", lambda=c(-3, 3), start=NULL, method="L-BFGS-B", ...) {
  data <- model.frame(object)
  if(any(unlist(lapply(as.list(data), class)) == "AsIs")) stop(
    "powerTransform for lmer models don't work with the 'I' function; rewrite your formula"
  )
  y <- (object@resp)$y
  fam <- match.fun(family)
  llik <- function(lambda){
    data$y.lambda <- fam(y, lambda, jacobian.adjusted=TRUE)
    m.lambda <- update(object, y.lambda ~ ., data=data)
    logLik(m.lambda)
  }
  if (is.null(start)) start <- 1
  res<- optimize(f = function(lambda1) llik(lambda1), lower=lambda[1], upper=lambda[2],
                 maximum=TRUE)
# optimize does not give the Hessian, so run optimHess
  res$hessian <- optimHess(res$maximum, llik, ...)
  res$invHess <- solve(-res$hessian)
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
  out <- data.frame(LRT=LR, df=df, pval=format.pval(pval))
  rownames(out) <-
    c(paste("LR test, lambda = (",
            paste(round(lambda, 2), collapse=" "), ")", sep=""))
  out}

