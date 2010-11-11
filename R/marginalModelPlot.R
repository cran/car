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
        main, AsIs=FALSE, ...){
  mf <- attr(model.frame(model), "terms")
  vform <- update(formula(model), terms)
  if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
     stop("Only predictors in the formula can be plotted. use mmp")
  terms <- attr(mf, "term.labels") # this is a list
  vterms <- attr(terms(vform), "term.labels")
# drop interactions (order > 1)
  vterms <- setdiff(vterms, terms[attr(mf, "order") > 1])
# keep only terms that are numeric or integer or factors or poly
  good <- NULL
  for (term in vterms) if(
      (AsIs == TRUE & inherits(model$model[[term]], "AsIs")) |
      inherits(model$model[[term]], "numeric") |
      inherits(model$model[[term]], "integer") |
      inherits(model$model[[term]], "poly")) good <- c(good, term)
  nt <- length(good) + fitted
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
  for (term in good){ 
    if(inherits(model$model[[term]], "poly")){
        horiz <- model.frame(model)[ , term][ , 1]
        lab <- paste("Linear part of", term)} else {
        horiz <- model.frame(model)[ , term]
        lab <- term }
    mmp(model, horiz, xlab=lab, ...)}
  if(fitted==TRUE) mmp(model, ...)
  mtext(side=3, outer=TRUE, main, line=0.5, cex=1.2)
  invisible()
  }
