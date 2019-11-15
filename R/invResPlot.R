# Last modified 25 Nov 2009 for point marking
# 18 January 2012 added robust estimation from Pendergast and Sheather
# 25 April 2016 check na.action for compatibility with Rcmdr
# 2017-02-13: modified to use id arg in calls to invTranPlot(). J. Fox
# 2017-11-30: substitute carPalette() for palette(). J. Fox
# 2019-05-16: make sure that xlab arg is properly passed to invTranPlot(). J. Fox
# 2019-11-14: change class(x) == "y" to inherits(x, "y")

inverseResponsePlot <- function(model, lambda=c(-1, 0, 1), robust=FALSE,
   xlab=NULL, ...)
    UseMethod("inverseResponsePlot")

invResPlot <- function(model, ...) UseMethod("inverseResponsePlot")


inverseResponsePlot.lm <- function(model, lambda=c(-1, 0, 1), robust=FALSE, 
       xlab=NULL, id=FALSE, ...) {
  if(inherits(model$na.action, "exclude")) model <- update(model, na.action=na.omit)
  id <- applyDefaults(id, defaults=list(method="x", n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
  if (isFALSE(id)){
      id.n <- 0
      id.method <- "none"
      labels <- id.cex <- id.col <- id.location <- NULL
  }
  else{
      labels <- id$labels
      if (is.null(labels)) labels <- names(residuals(model))
      id.method <- id$method
      id.n <- if ("identify" %in% id.method) Inf else id$n
      id.cex <- id$cex
      id.col <- id$col
      id.location <- id$location
  }
  if(robust == TRUE){
    m <- model$call
    m[[1L]] <- as.name("rlm")
    model <- eval(m, parent.frame())
  }
  mf <- model$model
  if (is.null(mf)) mf <- update(model, model=TRUE, method="model.frame")
  
  if (is.null(xlab)) xlab <- names(mf)[1] else force(xlab)
  y <- mf[, 1]
  yhat <- predict(model)
  invTranPlot(y, yhat, lambda=lambda, xlab=xlab, robust=robust, 
              id=list(n=id.n, method=id.method, labels=labels, cex=id.cex, col=id.col, location=id.location), ...)
}
