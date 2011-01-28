#############################################
# marginal model plots    Rev 12/30/09
# To do:
# Allow a Groups arg that will draw the plot for the specified group
# BUG:  sd's are WRONG with weights; see cards data
# 15 March 2010 changed to make
#   mmps(lm(longley)) work without specifying data or response
#   fixed bug  when only one plot is requested --- suppress call to par()
#   added 'outerLegend' to label lines
#   modified to work correctly with 
# 28 May 2010 S. Weisberg, fixed bugs in logistic models
#   changed line thickness of mean smooths
#   excluded SD smooth from bernoulli models
#   added grid lines
# 15 August 2010 fixed colors of points to work properly
# 16 January 2011 improved handling of splines and polynomials in mmps to
#    allow plots against base variables (e.g., bs(x, 3) could be
#    replaced by just x in the 'terms' argument to mmps.
#############################################
marginalModelPlot <- function(...){mmp(...)}
mmp <- function(model, ...){UseMethod("mmp")}

mmp.lm <- 
function (model, variable, mean = TRUE, sd = FALSE, 
    xlab = deparse(substitute(variable)), degree = 1, span = 2/3, key=TRUE,   
    ...)
{
    mmp.default(model, variable, mean, sd, xlab, degree, 
        span, key, ...)
}

mmp.default <-
function (model, variable, mean = TRUE, sd = FALSE, 
    xlab = deparse(substitute(variable)), degree = 1, span = 2/3, key=TRUE, 
    col.line = palette()[c(4, 2)], col=palette()[1], 
    labels, id.method="y", 
    id.n=if(id.method[1]=="identify") Inf else 0,
    id.cex=1, id.col=palette()[1], grid=TRUE, ...)
{   
    if  (!is.null(attr(model$model, "na.action"))) {
        if (attr(attr(model$model, "na.action"), "class") == "exclude")
            model <- update(model, na.action=na.omit)}     
    if (missing(variable)) {
        xlab <- "Fitted values"
        u <- fitted(model)
    } else {
        u <- variable}
    if(missing(labels)) 
        labels <- names(residuals(model))
    zpred <- function(...){pmax(predict(...), 0)}
    resp <- model.response(model.frame(model))
    plot(u, resp, 
        xlab = xlab, ylab = colnames(model$model[1]), type="n", ...)
	  if(grid){
      grid(lty=1, equilogs=FALSE)
      box()}
    points(u, model$model[ , 1], col=col, ...)
    if(key){
       outerLegend(c("Data", "Model"), lty=1:2, col=col.line, 
          bty="n", cex=0.75, fill=col.line, border=col.line, horiz=TRUE, 
          offset=0)
          }
    deg <- if(length(unique(u)) == 2) 0 else degree
    ow <- options(warn=-1)
    on.exit(options(ow))
    loess.y <- loess(resp ~ u, degree = deg, 
        span = span)
    loess.yhat <- loess(predict(model) ~ u, degree = deg, 
        span = span)
    new <- seq(min(u), max(u), length = 200)
    if (mean == TRUE) {
        lines(new, predict(loess.y, data.frame(u = new)), lty = 1, 
            lwd=2, col = col.line[1])
        lines(new, predict(loess.yhat, data.frame(u = new)), 
            lwd=2, lty = 2, col = col.line[2])
    }
    if (sd == TRUE) {
        loess.y.var <- loess(residuals(loess.y)^2 ~ u, degree = deg, 
            span = span)
        lines(new, predict(loess.y, data.frame(u = new)) +
            sqrt(zpred(loess.y.var, data.frame(u = new))), 
                 lty = 1, col = col.line[1])
        lines(new, predict(loess.y, 
            data.frame(u = new)) - sqrt(zpred(loess.y.var, 
            data.frame(u = new))), lty = 1, col = col.line[1])
        loess.yhat.var <- loess(residuals(loess.yhat)^2 ~ u, 
            degree = deg, span = span)
        s2 <- summary(model)$sigma^2
        lines(new, predict(loess.yhat, data.frame(u = new)) +
            sqrt(s2 + zpred(loess.yhat.var, data.frame(u = new))), 
            lty = 2, col = col.line[2])
        lines(new, predict(loess.yhat, data.frame(u = new)) -
            sqrt(s2 + zpred(loess.yhat.var, data.frame(u = new))), 
            lty = 2, col = col.line[2])
    }    
    showLabels(u, resp, labels=labels, 
        id.method=id.method, id.n=id.n, id.cex=id.cex, 
        id.col=id.col)
}

mmp.glm <- function (model, variable, mean = TRUE, sd = FALSE, 
    xlab = deparse(substitute(variable)), degree = 1, span = 2/3, key=TRUE, 
    col.line = palette()[c(4, 2)], col=palette()[1], 
    labels, id.method="y", 
    id.n=if(id.method[1]=="identify") Inf else 0,
    id.cex=1, id.col=palette()[1], grid=TRUE, ...)
{
    if (missing(variable)) {
        xlab <- "Linear Predictor"
        u <- fitted(update(model, na.action=na.omit))
    }  else {
        u <- variable }
    if(missing(labels)) 
        labels <- names(residuals(model)[!is.na(residuals(model))])
#    na.cases <- attr(model$model, "na.action")
#    if(length(na.cases)>0) u <- u[-na.cases]
    fr.mmp <- function(family, x) {
        if (family == "binomial")
            pmax(0, pmin(1, x))
        else if (family == "poisson")
            pmax(0, x)
        else if (family == "gamma")
            pmax(0, x)
        else x
    }
    response <- model.response(model.frame(model))
    fam <- model$family$family
    pw <- model$prior.weights # relevant only for binomial 
# For family = "binomial" we need to figure out the correct response
# The sample sizes are in  prior.weights
    bernoulli <- FALSE
    if(fam == "binomial") {
        if(!any(pw > 1.1)) bernoulli <- TRUE 
        if (is.factor(response)) {response <- as.numeric(response) - 1}
        if (is.matrix(response)){response <- response[, 1]/pw}  
    }
    plot(u, response, type="n", xlab = xlab, ylab = colnames(model$model[1]))
	  if(grid){
      grid(lty=1, equilogs=FALSE)
      box()}
    points(u, response, col=col, ...)
    if(key){
    outerLegend(c("Data", "Model"), lty=1:2, col=col.line, 
          bty="n", cex=0.75, fill=col.line, border=col.line, 
          horiz=TRUE, offset=0)
       }
    loess.y <- loess(response ~ u, degree = degree, span = span)
    loess.yhat <- loess(predict(model, type = "response") ~
        u, degree = degree, span = span)
    new <- seq(min(u), max(u), length = 200)
    pred.loess.y <- fr.mmp(fam, predict(loess.y, data.frame(u = new)))
    pred.loess.yhat <- fr.mmp(fam, predict(loess.yhat, data.frame(u = new)))
    if (mean == TRUE) {
        lines(new, pred.loess.y, lty = 1, col = col.line[1], lwd=2)
        lines(new, pred.loess.yhat, lty = 2, col = col.line[2], lwd=2)
    }
    if (sd == TRUE & bernoulli==FALSE) {
        loess.y.var <- loess(residuals(loess.y)^2 ~ u, degree = degree, 
            span = span)
        pred.loess.y.var <- pmax(0, predict(loess.y.var, data.frame(u = new)))
        lines(new, fr.mmp(fam, pred.loess.y + sqrt(pred.loess.y.var)), 
            lty = 1, col = col.line[1])
        lines(new, fr.mmp(fam, pred.loess.y - sqrt(pred.loess.y.var)), 
            lty = 1, col = col.line[1])
        loess.yhat.var <- loess(residuals(loess.yhat)^2 ~ u, 
            degree = degree, span = span)
        pred.loess.yhat.var <- pmax(0, predict(loess.yhat.var, 
            data.frame(u = new)))
        varfun <- summary(model)$dispersion * 
            model$family$variance(predict(model, type = "response"))/
              if (!is.null(model$prior.weights)) model$prior.weights else 1
        loess.varfun <- loess(varfun ~ u, degree = degree, span = span)
        pred.loess.varfun <- pmax(0, predict(loess.varfun, data.frame(u = new)))
        sd.smooth <- sqrt(pred.loess.yhat.var + pred.loess.varfun)
        lines(new, fr.mmp(fam, pred.loess.yhat + sd.smooth), 
            lty = 2, col = col.line[2])
        lines(new, fr.mmp(fam, pred.loess.yhat - sd.smooth), 
            lty = 2, col = col.line[2])
    }
    showLabels(u, model$model[, 1], labels=labels, 
        id.method=id.method, id.n=id.n, id.cex=id.cex, 
        id.col=id.col)
}

marginalModelPlots <- function(...) mmps(...)

mmps <- function(model, terms= ~ ., fitted=TRUE, layout=NULL, ask,
        main, ...){
  mf2 <- try(update(model, as.formula(terms), method="model.frame"),
     silent=TRUE)
# This second test is used for models like m1 <- lm(longley) which
# fail the first test becasue update doesn't work
  if(class(mf2) == "try-error")
       mf2 <- try(update(model, as.formula(terms),
               method="model.frame", data=model.frame(model)), silent=TRUE)
  if(class(mf2) == "try-error") stop("argument 'terms' not interpretable.")
  labels2 <- attr(attr(mf2, "terms"), "term.labels")
  order2 <- attr(attr(mf2, "terms"), "order")
  type2 <- rep("good", length(labels2))
  if(length(labels2) > 0) {
    for (j in 1:length(labels2)){
      if(order2[j] > 1) type2[j] <- NA #exclude interatctions
      if(inherits(mf2[[labels2[j]]], "factor")) type2[j] <- NA #no factors
      if(inherits(mf2[[labels2[j]]], "matrix")) type2[j] <- "original"
      }
    if (any( type2=="original", na.rm=TRUE )){
      p1 <- try(predict(model, type="terms"), silent=TRUE)
      if(class(p1) == "try-error") {type2[type2=="original"] <- NA} else
      warning("Splines and/or polynomials replaced by a fitted linear combination")
      }
  }
  nt <- sum(!is.na(type2)) + fitted
  if (missing(main)) main <- if (nt == 1) "Marginal Model Plot" else
     "Marginal Model Plots"
  if(is.null(layout)){
   layout <- switch(min(nt, 9),
           c(1, 1), c(1, 2), c(2, 2), c(2, 2), c(3, 2), c(3, 2),
           c(3, 3), c(3, 3), c(3, 3))}
  ask <- if(missing(ask) || is.null(ask)) prod(layout) < nt else ask
  if( prod(layout) > 1) {
    op <- par(mfrow=layout, ask=ask, no.readonly=TRUE,
            oma=c(0, 0, 2.5, 0), mar=c(5, 4, 1.5, 1.5) + .1)
    on.exit(par(op))
  }
  if (length(labels2) > 0) {
    for (j in 1:length(labels2)) {
      if(!is.na(type2[j])) {
        horiz <- if(type2[j] == "original"){p1[, labels2[j]]} else {
                 if(type2[j] == "good") mf2[ , labels2[j]] else NULL}
        lab <- labels2[j]
        mmp(model, horiz, xlab=lab, ...)}}
    }
  if(fitted==TRUE) mmp(model, ...)
  if(nt==1)
     mtext(side=3, outer=FALSE, main, line=1.5, cex=1.2) else
     mtext(side=3, outer=TRUE,  main, line=0.5, cex=1.2)
  if(any(is.na(type2))) warning("Interactions and/or factors skipped")
  invisible()
  }
