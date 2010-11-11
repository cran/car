# Modified Nov. 24, 2009 by S. Weisberg to use showLabels
#   rather than showExtremes
# 11 & 20 January 2010: changed lty=3 to lty=1 for fitted curve. J. Fox
# 14 April 2010: set id.n = 0. J. Fox
# 15 April 2010; rewrite showLabels
# 25 May 2010 added grid() to plots, S. Weisberg
# 15 August 2010, fixed so col= works correctly with plot, but not Boxplot
# 15 August 2010, deleted pch= argument, as it wasn't used

residualPlots <- function(model, ...){UseMethod("residualPlots")}

residualPlots.default <- function(model, terms= ~ . , 
     layout=NULL, ask, main="", 
     fitted=TRUE, AsIs=FALSE, plot=TRUE, tests=TRUE, ...){
  mf <- attr(model.frame(model), "terms")
  vform <- update(formula(model), terms)
  if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
     stop("Only predictors in the formula can be plotted.")
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
      inherits(model$model[[term]], "factor") | 
      inherits(model$model[[term]], "poly")) good <- c(good, term)
  nt <- length(good) + fitted
  if (nt == 0) stop("No plots specified")
  if(is.null(layout)){
   layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 2), c(2, 2), 
                               c(3, 2), c(3, 2), c(3, 3), c(3, 3), c(3, 3))}
  nr <- 0
  ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
  if(nt > 1 & plot == TRUE) {
    op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
            oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
    on.exit(par(op))}
  ans <- NULL       
  if(!is.null(good)){
    for (term in good){
     nr <- nr + 1
     qtest <- residualPlot(model, term, plot=plot, ...)
     if(!is.null(qtest)){
        ans <- rbind(ans, qtest)
        row.names(ans)[nr] <- term}
    } }   
  # Tukey's test
  if (fitted == TRUE){      
   tuk <- residualPlot(model, "fitted", plot=plot, ...)
   if (!is.null(tuk)  & class(model)[1] == "lm"){
      ans <- rbind(ans, tuk)
      row.names(ans)[nr + 1] <- "Tukey test"
      ans[nr + 1, 2] <- 2*pnorm(abs(ans[nr + 1, 1]), lower.tail=FALSE)}} 
  if(plot == TRUE) mtext(side=3, outer=TRUE, main, cex=1.2)
  if(!is.null(ans)) {
     dimnames(ans)[[2]] <- c("Test stat", "Pr(>|t|)")
     return(if(tests == FALSE) invisible(ans) else round(ans, 3)) } else
  invisible(NULL)
  }
  
residualPlots.lm <- function(model, ...) {
 residualPlots.default(model, ...)
 }
  
residualPlots.glm <- function(model, ...) {
 residualPlots.default(model, ...)
 }

residualPlot <- function(model, ...) UseMethod("residualPlot")

residualPlot.default <- function(model, variable = "fitted", type = "pearson", 
                 plot = TRUE,     
                 quadratic = TRUE, 
                 smooth = FALSE, span = 1/2, smooth.lwd=lwd, smooth.lty=lty, 
                 smooth.col=col.lines,  
                 labels, 
                 id.method = "xy", 
                 id.n = if(id.method[1]=="identify") Inf else 0,
                 id.cex=1, id.col=palette()[1], 
                 col = palette()[1], col.lines = palette()[2], 
                 xlab, ylab, lwd = 1, lty = 1,  
                 grid=TRUE, ...) {
# two functions modified from 'scatterplot' function:
	logged <- function(axis=c("x", "y")){
		axis <- match.arg(axis)
		0 != length(grep(axis, log))
	}
	lowess.line <- function(x, y, col, span) {
	  call <- match.call(expand.dots=TRUE)
	  uselogx <- length(call$log == "x") == 1
    if(uselogx) x <- log(x)
		valid <- complete.cases(x, y)
		x <- x[valid]
		y <- y[valid]
		ord <- order(x)
		x <- x[ord]
		y <- y[ord]
		fit <- loess.smooth(x, y, span=span)
		lines(if(uselogx) exp(fit$x) else fit$x, fit$y, 
         lwd=smooth.lwd, col=smooth.col, lty=smooth.lty)
		}
# End of copied functions
 string.capitalize <- function(string) {
     paste(toupper(substring(string, 1, 1)), substring(string, 2), sep="")}
 if(missing(labels)) 
      labels <-  names(residuals(model)[!is.na(residuals(model))])
 ylab <- if(!missing(ylab)) ylab else
         paste(string.capitalize(type), "residuals")
 column <- match(variable, names(model$model))
 if(is.na(column) && variable != "fitted")
   stop(paste(variable, "is not a term in the mean function"))
 horiz <- if(variable == "fitted") predict(model) else model$model[[column]]
 lab <- if(variable == "fitted") {
    if(inherits(model, "glm")) 
       "Linear Predictor" else "Fitted values"} else variable
 lab <- if(!missing(xlab)) xlab else lab
 ans <-
   if(inherits(horiz, "poly")) {
       horiz <- horiz[ , 1]
       lab <- paste("Linear part of", lab)
       c(NA, NA)}
   else if (class(horiz) == "factor") c(NA, NA)
   else residCurvTest(model, variable)
# ans <- if (class(horiz) != "factor")  else c(NA, NA)
 if(plot==TRUE){
  vert <- switch(type, "rstudent"=rstudent(model), 
       "rstandard"=rstandard(model), residuals(model, type=type))
  if(class(horiz) == "factor") {
     idm <- switch(id.method, xy="y", x="y", y="y", "none")  
     Boxplot(vert, horiz, xlab=lab, ylab=ylab, labels=labels, 
            id.method=idm, id.n=id.n, id.cex=id.cex,  
            id.col=id.col, ...) 
     abline(h=0, lty=2) } else {    
     plot(horiz, vert, xlab=lab, ylab=ylab, ...)
	  if(grid){
      grid(lty=1, equilogs=FALSE)
      box()}
     points(horiz, vert, col=col, ...)
     abline(h=0, lty=2)
     if(quadratic==TRUE){
        new <- seq(min(horiz), max(horiz), length=200)
        lm2 <- lm(residuals(model, type="pearson")~poly(horiz, 2))
        lines(new, predict(lm2, list(horiz=new)), lty=1, lwd=2, col=col.lines)
        }
     if(smooth==TRUE) lowess.line(horiz, vert, smooth.col, span)
     showLabels(horiz, vert, labels=labels, 
            id.method=id.method, id.n=id.n, id.cex=id.cex, 
            id.col=id.col)  
        }
      }  
  invisible(ans)}
 
residCurvTest <- function(model, variable) {UseMethod("residCurvTest")}
residCurvTest.lm <- function(model, variable) {
 if(variable == "fitted") tukeyNonaddTest(model) else {
  if(is.na(match(variable, attr(model$terms, "term.labels"))))
     stop(paste(variable, "is not a term in the mean function")) else {
     xsqres <- qr.resid(model$qr, model.frame(model)[[variable]]^2)
     r <- residuals(model, type="pearson")
     m1 <- lm(r ~ xsqres, weights=weights(model))
     df.correction <- sqrt((df.residual(model)-1) / df.residual(m1))
     test <- summary(m1)$coef[2, 3] * df.correction
     c(Test=test, Pvalue=2 * pt(-abs(test), df.residual(model)-1))
     }}}

residCurvTest.glm <- function(model, variable) {
 if(variable == "fitted") c(NA, NA) else {
  if(is.na(match(variable, attr(model$terms, "term.labels"))))
     stop(paste(variable, "is not a term in the mean function")) else {
     newmod <- paste(" ~ . + I(", variable, "^2)", sep="")
     m2 <- update(model, newmod)
     c(Test= test<-deviance(model)-deviance(m2), Pvalue=1-pchisq(test, 1))
}}}
     
tukeyNonaddTest <- function(model){
 qr <- model$qr
 fitsq <- qr.resid(qr, predict(model, type="response")^2)
 r <- residuals(model, type="pearson")
 m1 <- lm(r~fitsq, weights=weights(model))
 df.correction <- sqrt((df.residual(model) - 1)/df.residual(m1))
 tukey <- summary(m1)$coef[2, 3] * df.correction
 c(Test=tukey, Pvalue=2*pnorm(-abs(tukey)))
 }
 
residualPlot.lm <- function(model, ...) {
  residualPlot.default(model, ...)
  }
  
residualPlot.glm <- function(model, variable = "fitted", type = "pearson", 
                 plot = TRUE, quadratic = FALSE, smooth = TRUE, ...) {
  residualPlot.default(model, variable=variable, type=type, plot=plot, 
                 quadratic=quadratic, smooth=smooth, ...)
  }

