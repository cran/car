# October 23, 2009  avPlots by S. Weisberg.  avPlot by John Fox
# 13 January 2010: changed default id.n=3. J. Fox
# 13 March 2010: added intercept argument. J. Fox
# 14 April 2010: set id.n = 0. J. Fox
# 22 April 2010: modified id.n S. Weisberg
# 10 May 2010:  added gridlines
# 25 May 2010:  changed default color scheme

avPlots <- function(model, terms=~., intercept=FALSE, layout=NULL, ask, 
           main, ...){
  terms <- if(is.character(terms)) paste("~",terms) else terms
  vform <- update(formula(model),terms)
  if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
     stop("Only predictors in the formula can be plotted.")
  terms.model <- attr(attr(model.frame(model), "terms"), "term.labels")
  terms.vform <- attr(terms(vform), "term.labels")
  terms.used <- match(terms.vform, terms.model)
  mm <- model.matrix(model) 
  model.names <- attributes(mm)$dimnames[[2]]
  model.assign <- attributes(mm)$assign
  good <- model.names[!is.na(match(model.assign, terms.used))]
  if (intercept) good <- c("(Intercept)", good)
  nt <- length(good)
  if (nt == 0) stop("No plots specified")
  if (missing(main)) main <- if (nt == 1) "Added-Variable Plot" else "Added-Variable Plots"
  if(is.null(layout)){
   layout <- switch(min(nt,9), c(1,1), c(1,2), c(2,2), c(2,2),
                    c(3,2), c(3,2), c(3,3), c(3,3), c(3,3))}
  nr <- 0
  ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
  if(nt > 1) {
    op<-par(no.readonly=TRUE, oma=c(0, 0, 1.5, 0), 
          mar=c(5, 4, 1, 2) + .1, mfrow=layout, ask=ask)
    on.exit(par(op))}
  for (term in good) avPlot(model, term, main="", ...)
  mtext(side=3,outer=TRUE,main, cex=1.2)
  invisible(NULL)
 }
 
avp <- function(...) avPlots(...)
 
avPlot <-  function(model, ...) UseMethod("avPlot")

avPlot.lm <-
function (model, variable,
    id.method = list(abs(residuals(model, type="pearson")), "x"),
    labels, 
    id.n = if(id.method[1]=="identify") Inf else 0,
    id.cex=1, id.col=palette()[1],
    col = palette()[1], col.lines = palette()[2],
    xlab, ylab, pch = 1, lwd = 2, main="Added-variable Plot", grid=TRUE, ...)
{
    variable <- if (is.character(variable) & 1 == length(variable))
        variable
    if(missing(labels)) 
        labels <- names(residuals(model)[!is.na(residuals(model))])
    else deparse(substitute(variable))
    mod.mat <- model.matrix(model)
    var.names <- colnames(mod.mat)
    var <- which(variable == var.names)
    if (0 == length(var))
        stop(paste(variable, "is not a column of the model matrix."))
    response <- response(model)
    responseName <- responseName(model)
    if (is.null(weights(model)))
        wt <- rep(1, length(response))
    else wt <- weights(model)
    res <- lsfit(mod.mat[, -var], cbind(mod.mat[, var], response),
        wt = wt, intercept = FALSE)$residuals
    xlab <- if(missing(xlab)) paste(var.names[var], "| others") else xlab
    ylab <- if(missing(ylab)) paste(responseName, " | others")  else ylab
    plot(res[, 1], res[, 2], xlab = xlab, ylab = ylab, type="n", ...)
	  if(grid){
      grid(lty=1, equilogs=FALSE)
      box()}
    points(res[, 1], res[, 2], col=col, pch=pch, ...)
    abline(lsfit(res[, 1], res[, 2], wt = wt), col = col.lines, lwd = lwd)
    showLabels(res[, 1],res[, 2], labels=labels, 
          id.method=id.method, id.n=id.n, id.cex=id.cex, 
          id.col=id.col)  
}

avPlot.glm<-function(model, variable, 
    id.method = list(abs(residuals(model, type="pearson")), "x"),
    labels,
    id.n = if(id.method[1]=="identify") Inf else 0,
    id.cex=1, id.col=palette()[1], 
    col = palette()[1], col.lines = palette()[2],
    xlab, ylab, pch = 1, lwd = 2,  type=c("Wang", "Weisberg"), 
    main="Added-variable Plot", grid=TRUE, ...){
    #last modified 20 Feb 2002 by J. Fox
    type<-match.arg(type)
    if(missing(labels)) labels <- names(residuals(model)[!is.na(residuals(model))])
    variable<-if (is.character(variable) & 1==length(variable)) variable
        else deparse(substitute(variable))
    mod.mat<-model.matrix(model)
    var.names<-colnames(mod.mat)
    var<-which(variable==var.names)
    if (0==length(var)) stop(paste(variable,"is not a column of the model matrix."))
    response<-response(model)
    responseName<-responseName(model)
    wt<-model$prior.weights
    mod<-glm(response~mod.mat[,-var]-1, weights=wt, family=family(model))
    res.y<-residuals(mod, type="pearson")
    wt<-if (type=="Wang") wt*model$weights else wt
    res.x<-lsfit(mod.mat[,-var], mod.mat[,var], wt=wt,    
        intercept=FALSE)$residuals
    xlab <- if(missing(xlab)) paste(var.names[var], "| others") else xlab
    ylab <- if(missing(ylab)) paste(responseName, " | others") else ylab
    plot(res.x, res.y, xlab=xlab, type="n", ylab=ylab, main=main, ...)
	  if(grid){
      grid(lty=1, equilogs=FALSE)
      box()}
    points(res.x, res.y, col=col, pch=pch, ...)
    abline(lsfit(res.x, res.y, wt=wt), col=col.lines, lwd=lwd)
    showLabels(res.x,res.y, labels=labels, 
          id.method=id.method, id.n=id.n, id.cex=id.cex, 
          id.col=id.col)  
    }

  
  
