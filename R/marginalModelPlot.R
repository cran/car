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
# 16 June 2011 allow layout=NA, in which case the layout is not set in this
#  function, so it is the responsibility of the user
# 14 September 2012 improved the smoothers
#############################################
marginalModelPlot <- function(...){mmp(...)}
mmp <- function(model, ...){UseMethod("mmp")}

mmp.lm <- 
function (model, variable, sd = FALSE,
    xlab = deparse(substitute(variable)),
    smoother = loessLine, smoother.args=list(degree=1, span=2/3),
    key=TRUE, ...)
{
    mmp.default(model, variable, sd, xlab, smoother=smoother,
       smoother.args=smoother.args, key, ...)
}

mmp.default <-
function (model, variable, sd = FALSE,
    xlab = deparse(substitute(variable)), smoother=loessLine,
    smoother.args, key=TRUE,
    col.line = palette()[c(4, 2)], col=palette()[1], 
    labels, id.method="y", 
    id.n=if(id.method[1]=="identify") Inf else 0,
    id.cex=1, id.col=palette()[1], grid=TRUE, ...)
{     
    if (!is.function(smoother)) {
        smoother <- loessLine
        smoother.args <- list()
        }
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
    ow <- options(warn=-1)
    on.exit(options(ow))
    smoother.args$lty <- smoother.args$lty.spread <- 1
    smoother(u, resp, col.line[1], log.x=FALSE, log.y=FALSE, spread=sd,
       smoother.args=smoother.args)
    smoother.args$lty <- smoother.args$lty.spread <- 2
    smoother(u, predict(model), col.line[2], log.x=FALSE, log.y=FALSE, spread=sd,
       smoother.args=smoother.args)
    showLabels(u, resp, labels=labels, 
        id.method=id.method, id.n=id.n, id.cex=id.cex, 
        id.col=id.col)
}

mmp.glm <- function (model, variable, sd = FALSE, 
    xlab = deparse(substitute(variable)), smoother=gamLine,
    smoother.args=list(k=3), key=TRUE, 
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
    response <- model.response(model.frame(model))
    fam <- model$family$family
    lin <- model$family$link
    pw <- model$prior.weights # relevant only for binomial 
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
    smoother.args$lty <- smoother.args$lty.spread <- 1
    smoother.args$family <- fam
    smoother.args$link <- lin
    smoother.args$weights <- pw
    smoother(u, response, col.line[1], log.x=FALSE, log.y=FALSE, spread=FALSE,
       smoother.args=smoother.args)
    smoother.args$lty <- smoother.args$lty.spread <- 2
    response <- if(fam=="binomial") predict(model, type="response")/pw
         else predict(model, type="response")
    smoother(u, response, col.line[2], log.x=FALSE, log.y=FALSE, spread=FALSE,
       smoother.args=smoother.args)
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
# fail the first test because update doesn't work
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
  if (nt == 0) stop("No plots specified")
  if (nt > 1 & (is.null(layout) || is.numeric(layout))) {
    if(is.null(layout)){
         layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 2), c(2, 2), 
                             c(3, 2), c(3, 2), c(3, 3), c(3, 3), c(3, 3))
    }
    ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
    op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
            oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
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
