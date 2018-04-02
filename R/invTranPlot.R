# Modified 25 Nov 2009 for point marking
# 20 Jan 2010: changed line types. J. Fox
# 15 August 2010: fixed colors of points
# 18 January 2011; added robust M estimation
# 2017-02-13: consolidated id argument. John
# 2017-02-16: replace rlm() calls with MASS::rlm()
# 2017-11-30: substitute carPalette() for palette(). J. Fox

invTranPlot <- function(x,...) UseMethod("invTranPlot")

invTranPlot.formula <- function(x, data, subset, na.action, id=FALSE, ...) {
    mf <- call <- match.call()
    m <- match(c("x", "data", "subset", "na.action"), 
         names(mf), 0L)
    mf <- mf[c(1L,m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")     
    names(mf)[which(names(mf)=="x")] <- "formula"
    mf <- eval(mf, parent.frame())
    if(dim(mf)[2] != 2) stop("Formula must be of the form y ~ x")
    id <- applyDefaults(id, defaults=list(method="x", n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
    if (isFALSE(id)){
        id.n <- 0
        id.method <- "none"
        labels <- id.cex <- id.col <- id.location <- NULL
    }
    else{
        labels <- id$labels
        if (is.null(labels)) labels <- rownames(mf)
        id.method <- id$method
        id.n <- if ("identify" %in% id.method) Inf else id$n
        id.cex <- id$cex
        id.col <- id$col
        id.location <- id$location
    }
    vx <- mf[,2]
    vy <- mf[,1]
    if( is.null(call$xlab) & is.null(call$ylab)) 
         invTranPlot(vx, vy, xlab=colnames(mf)[2], ylab=colnames(mf)[1], 
                     id=list(n=id.n, method=id.method, labels=labels, cex=id.cex, col=id.col, location=id.location), ...) else
     if(is.null(call$xlab) & !is.null(call$ylab)) 
         invTranPlot(vx, vy, xlab=colnames(mf)[2], 
                     id=list(n=id.n, method=id.method, labels=labels, cex=id.cex, col=id.col, location=id.location), ...) else
     if(!is.null(call$xlab) & is.null(call$ylab)) 
         invTranPlot(vx ,vy ,ylab=colnames(mf)[1], 
                     id=list(n=id.n, method=id.method, labels=labels, cex=id.cex, col=id.col, location=id.location), ...) else
         invTranPlot(vx, vy, id=list(n=id.n, method=id.method, labels=labels, cex=id.cex, col=id.col, location=id.location), ...)
  }   

invTranPlot.default<- function(x, y, lambda=c(-1, 0, 1), robust=FALSE, 
        lty.lines=rep(c("solid", "dashed", "dotdash", "longdash", "twodash"), 
        length=1 + length(lambda)), lwd.lines=2, 
        col=carPalette()[1], col.lines=carPalette(), 
        xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
        family="bcPower", optimal=TRUE, key="auto", id=FALSE,
        grid=TRUE, ...){
    id <- applyDefaults(id, defaults=list(method="x", n=2, cex=1, col=carPalette()[1], location="lr"), type="id")
    if (isFALSE(id)){
        id.n <- 0
        id.method <- "none"
        labels <- id.cex <- id.col <- id.location <- NULL
    }
    else{
        labels <- id$labels
        if (is.null(labels)) labels <- seq(along=x)
        id.method <- id$method
        id.n <- if ("identify" %in% id.method) Inf else id$n
        id.cex <- id$cex
        id.col <- id$col
        id.location <- id$location
    }
 if (is.factor(x)) stop("Predictor variable may not be a factor")
 if (is.factor(y)) stop("Response variable may not be a factor")
 if (optimal){
     opt <- invTranEstimate(x, y, family=family, confidence=FALSE, robust=robust)
     lam <- c(opt$lambda, lambda)} else lam <- lambda
 fam <- match.fun(family)
 plot(x, y, xlab=xlab, ylab=ylab, type="n", col=col, ...)
 if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
 points(x, y, col=col, ...)
 rss <- NULL
 new <- seq(min(x, na.rm=TRUE), max(x,na.rm=TRUE), length=100)
 for (j in 1:length(lam)){
     m1 <- if(robust) MASS::rlm(y ~ fam(x, lam[j])) else lm(y~fam(x, lam[j]))
     rss <- c(rss, sum(residuals(m1)^2))
     lines(new,predict(m1, data.frame(x=new)), lty=lty.lines[j],
       col=col.lines[j], lwd=lwd.lines)}
 showLabels(x, y, labels=labels, 
          method=id.method, n=id.n, cex=id.cex, 
          col=id.col, location = id.location)
 if (!is.null(key)) {
      loc <- key
      if(length(lam) <= 4) {
        lims <- par("usr")[c(1,4)]
        llam <- expression(paste(hat(lambda), ":"))
        text(lims[1],lims[2], llam, xpd=TRUE, pos=3)
        outerLegend(
            as.character(round(lam,2)),
            lwd=lwd.lines, lty=lty.lines, 
            col=col.lines,
            bty="n", 
            cex=0.85, fill=col.lines, 
            border=col.lines, horiz=TRUE, adjust=FALSE)}
      else {
        legend(ifelse(cor(x, y)>0,"bottomright","topright"), 
            legend =  c(expression(hat(lambda)),as.character(round(lam,2))), 
            lwd=lwd.lines, lty=c("blank", lty.lines), col=c("#00000000",col.lines), 
            inset=0.02,  cex=0.75, fill=c("#00000000",col.lines), 
            border=c("#00000000",col.lines))
      }}
 data.frame(lambda=lam, RSS=rss)
}

invTranEstimate <- function(x, y, family="bcPower", confidence=0.95,
  robust=FALSE){
  if (is.factor(x)) stop("Predictor variable may not be a factor")
  if (is.factor(y)) stop("Response variable may not be a factor")
  if (robust) confidence <- FALSE
  fam <- match.fun(family)
  f <- if(robust==FALSE)
    function(lambda,x,y,family){deviance(lm(y~fam(x,lambda)))}  else
    function(lambda,x,y,family){sum(residuals(MASS::rlm(y ~ fam(x,lambda)))^2)}
  lhat <- optimize(f = function(lambda) f(lambda, x, y, family),interval=c(-10,10))
  if (confidence==FALSE){ return(list(lambda=lhat$minimum)) } else {
    g <- lm(y~fam(x,lhat$minimum))
    n = length(residuals(g))
    dev0 <- -n*log(deviance(g))
    cutoff <- qchisq(confidence,1)/2
    f1 <- function(lam) abs(dev0 + n*log(deviance(lm(y~fam(x,lam)))) -cutoff)
    lowlim <- optimize(f1, interval=c(-10,lhat$minimum))
    hilim <-  optimize(f1, interval=c(lhat$minimum,10))
    return(list(lambda=lhat$minimum,lowerCI=lowlim$minimum,upperCI=hilim$minimum))}
}
        

