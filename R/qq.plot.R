# Quantile-comparison plots (J. Fox)

# last modified 8 October 05 by J. Fox

qqp<-function(...) qq.plot(...)

qq.plot<-function(x, ...) {
    UseMethod("qq.plot")
    }
  
qq.plot.default<-function(x, distribution="norm", ylab=deparse(substitute(x)),
        xlab=paste(distribution, "quantiles"), main=NULL, las=par("las"),
        envelope=.95, labels=FALSE, col=palette()[2], lwd=2, pch=1, cex=1,
        line=c("quartiles", "robust", "none"), ...){
    # last modified 23 February 2003
    result <- NULL
    line<-match.arg(line)
    good<-!is.na(x)
    ord<-order(x[good])
    ord.x<-x[good][ord]
    q.function<-eval(parse(text=paste("q",distribution, sep="")))
    d.function<-eval(parse(text=paste("d",distribution, sep="")))
    n<-length(ord.x)
    P<-ppoints(n)
    z<-q.function(P, ...)
    plot(z, ord.x, xlab=xlab, ylab=ylab, main=main, las=las, col=col, pch=pch,
        cex=cex)
    if (line=="quartiles"){
        Q.x<-quantile(ord.x, c(.25,.75))
        Q.z<-q.function(c(.25,.75), ...)
        b<-(Q.x[2]-Q.x[1])/(Q.z[2]-Q.z[1])
        a<-Q.x[1]-b*Q.z[1]
        abline(a, b, col=col, lwd=lwd)
        }
    if (line=="robust"){
        if (!require("MASS")) stop("MASS package not available")
        coef<-coefficients(rlm(ord.x~z))
        a<-coef[1]
        b<-coef[2]
        abline(a,b)
        }
    if (line != 'none' & envelope != FALSE) {
        zz<-qnorm(1-(1-envelope)/2)
        SE<-(b/d.function(z, ...))*sqrt(P*(1-P)/n)
        fit.value<-a+b*z
        upper<-fit.value+zz*SE
        lower<-fit.value-zz*SE
        lines(z, upper, lty=2, lwd=lwd/2, col=col)
        lines(z, lower, lty=2, lwd=lwd/2, col=col)
        }
    if (labels[1]==TRUE & length(labels)==1) labels<-seq(along=z)
    if (labels[1] != FALSE) {
        selected<-identify(z, ord.x, labels[good][ord])
        result <- seq(along=x)[good][ord][selected]
        }
    if (is.null(result)) invisible(result) else sort(result)
    }
    
qq.plot.lm<-function(x, main=NULL, xlab=paste(distribution, "Quantiles"),
    ylab=paste("Studentized Residuals(",deparse(substitute(x)),")",sep=""),
    distribution=c("t", "norm"), line=c("quartiles", "robust", "none"), las=par("las"),
    simulate=FALSE, envelope=.95, labels=names(rstudent), reps=100, 
    col=palette()[2], lwd=2, pch=1, cex=1, ...){
    # last modified 23 Feb 2003
    result <- NULL
    distribution <- match.arg(distribution)
    line<-match.arg(line)
    rstudent<-rstudent(x)
    sumry <- summary.lm(x)
    res.df<-sumry$df[2]
    if(!simulate){
        if (distribution == 't')
            result <- qq.plot.default(rstudent, distribution='t', df=res.df-1, line=line,
                main=main, xlab=xlab, ylab=ylab, las=las, envelope=envelope, labels=labels, 
                col=col, lwd=lwd, pch=pch, cex=cex, ...)
        else
            result <- qq.plot.default(rstudent, distribution='norm', line=line,
                main=main, xlab=xlab, ylab=ylab, las=las, envelope=envelope, labels=labels, 
                col=col, lwd=lwd, pch=pch, cex=cex, ...) 
        }
    else {
        good <- !is.na(rstudent)
        n<-length(rstudent)
        rstudent <- na.omit(rstudent)
        ord<-order(rstudent)
        ord.x<-rstudent[ord]
        n<-length(ord)
        P<-ppoints(n)
        z<-if (distribution == 't') qt(P, df=res.df-1) else qnorm(P)
        plot(z, ord.x, xlab=xlab, ylab=ylab, main=main, las=las, pch=pch, col=col, cex=cex)
        yhat<-na.omit(fitted.values(x))
        S<-sumry$sigma
        Y<-matrix(yhat,n,reps)+matrix(rnorm(n*reps, sd=S),n,reps)
        X<-model.matrix(x)
        rstud<-apply(rstudent(lm(Y~X-1)),2,sort)
        lower<-apply(rstud,1,quantile,prob=(1-envelope)/2)
        upper<-apply(rstud,1,quantile,prob=(1+envelope)/2)
        lines(z, upper, lty=2, lwd=lwd/2, col=col)
        lines(z, lower, lty=2, lwd=lwd/2, col=col)
        if (line=="quartiles"){
            Q.x<-quantile(rstudent, c(.25,.75))
            Q.z <- if (distribution == 't') qt(c(.25,.75),df=res.df-1) else qnorm(c(.25,.75))
            b<-(Q.x[2]-Q.x[1])/(Q.z[2]-Q.z[1])
            a<-Q.x[1]-b*Q.z[1]
            abline(a, b, col=col, lwd=lwd)
            }
        if (line=="robust"){
            if (!require("MASS")) stop("MASS package not available")
            coef<-coefficients(rlm(ord.x~z))
            a<-coef[1]
            b<-coef[2]
            abline(a, b, col=col, lwd=lwd)
            }
        if (labels[1]==TRUE & length(labels)==1) labels<-seq(along=z)
        if (labels[1] != FALSE) {
            selected<-identify(z, ord.x, labels[ord])
            result <- (1:n)[good][ord][selected]
            }
        }
    if (is.null(result)) invisible(result) else sort(result)
    }
 
qq.plot.glm<-function(x, ...){
    stop("QQ plot for studentized residuals not available for glm")
    }
