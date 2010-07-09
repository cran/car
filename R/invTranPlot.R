# Modified 25 Nov 2009 for point marking
# 20 Jan 2010: changed line types. J. Fox

invTranPlot <- function(x,...) UseMethod("invTranPlot")

invTranPlot.formula <- function(x, data, subset, na.action, ...) {
    mf <- call <- match.call()
    m <- match(c("x", "data", "subset", "na.action"), 
         names(mf), 0L)
    mf <- mf[c(1L,m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")     
    names(mf)[which(names(mf)=="x")] <- "formula"
    mf <- eval(mf, parent.frame())
    if(dim(mf)[2] != 2) stop("Formula must be of the form y ~ x")
    vx <- mf[,2]
    vy <- mf[,1]
    if( is.null(call$xlab) & is.null(call$ylab)) 
         invTranPlot(vx,vy,xlab=colnames(mf)[2],ylab=colnames(mf)[1],...) else
     if(is.null(call$xlab) & !is.null(call$ylab)) 
         invTranPlot(vx,vy,xlab=colnames(mf)[2],...) else
     if(!is.null(call$xlab) & is.null(call$ylab)) 
         invTranPlot(vx,vy,ylab=colnames(mf)[1],...) else
         invTranPlot(vx,vy,...)
  }   

invTranPlot.default<- function(x, y, lambda=c(-1, 0, 1), 
        lty.lines=rep(c("solid", "dashed", "dotdash", "longdash", "twodash"), 
        length=1 + length(lambda)), lwd.lines=2, 
        col.lines=palette(),
        xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
        family="bcPower", optimal=TRUE, key="auto",
        id.method = abs(residuals(lm(y~x))),
        labels, 
        id.n = if(id.method[1]=="identify") Inf else 0, 
        id.cex=1, id.col=palette()[1], grid=TRUE, ...){
 if (missing(labels)) labels <- seq(length(x))
 if (is.factor(x)) stop("Predictor variable may not be a factor")
 if (is.factor(y)) stop("Response variable may not be a factor")
 if (optimal){opt <- invTranEstimate(x, y, family=family, confidence=FALSE)
              lam <- c(opt$lambda, lambda)} else lam <- lambda
 fam <- match.fun(family)
 plot(x, y, xlab=xlab, ylab=ylab, type="n", ...)
 if(grid){
    grid(lty=1, equilogs=FALSE)
    box()}
 points(x, y, ...)
 rss <- NULL
 new <- seq(min(x, na.rm=TRUE), max(x,na.rm=TRUE), length=100)
 for (j in 1:length(lam)){
     m1 <- lm(y~fam(x, lam[j]))
     rss <- c(rss, deviance(m1))
     lines(new,predict(m1, data.frame(x=new)), lty=lty.lines[j],
       col=col.lines[j], lwd=lwd.lines)}
 showLabels(x, y, labels=labels, 
          id.method=id.method, id.n=id.n, id.cex=id.cex, 
          id.col=id.col)
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
 data.frame(lambda=lam,RSS=rss)
}

invTranEstimate <- function(x,y,family="bcPower",confidence=0.95){
  if (is.factor(x)) stop("Predictor variable may not be a factor")
  if (is.factor(y)) stop("Response variable may not be a factor")
  fam <- match.fun(family)
  f <- function(lambda,x,y,family){deviance(lm(y~fam(x,lambda)))}
  lhat <- optimize(f = function(lambda) f(lambda,x,y,family),interval=c(-10,10))
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
        

