# Quantile-comparison plots (J. Fox)

qqp<-function(...) qq.plot(...)

qq.plot<-function(x, ...) {
    UseMethod("qq.plot")
    }
  
qq.plot.default<-function(x, distribution="norm", ylab=deparse(substitute(x)),
        xlab=paste(distribution, "quantiles"), main="", las=1,
        envelope=.95, labels=F, col=palette()[2], lwd=2, pch=1,
        line=c("quartiles", "robust"), ...){
    # last modified 1 Feb 2001
    line<-match.arg(line)
    good<-!is.na(x)
    ord<-order(x[good])
    ord.x<-x[good][ord]
    q.function<-eval(parse(text=paste("q",distribution, sep="")))
    d.function<-eval(parse(text=paste("d",distribution, sep="")))
    n<-length(ord.x)
    P<-ppoints(n)
    z<-q.function(P, ...)
    plot(z, ord.x, xlab=xlab, ylab=ylab, main=main, las=las, col=col, pch=pch)
    if (line=="quartiles"){
        Q.x<-quantile(ord.x, c(.25,.75))
        Q.z<-q.function(c(.25,.75), ...)
        b<-(Q.x[2]-Q.x[1])/(Q.z[2]-Q.z[1])
        a<-Q.x[1]-b*Q.z[1]
        abline(a, b, col=col, lwd=lwd)
        }
    if (line=="robust"){
        if (!require("MASS", quietly=T)) stop("MASS package not available")
        coef<-coefficients(rlm(ord.x~z))
        a<-coef[1]
        b<-coef[2]
        abline(a,b)
        }
    if (envelope != F) {
        zz<-qnorm(1-(1-envelope)/2)
        SE<-(b/d.function(z, ...))*sqrt(P*(1-P)/n)
        fit.value<-a+b*z
        upper<-fit.value+zz*SE
        lower<-fit.value-zz*SE
        lines(z, upper, lty=2, lwd=lwd/2, col=col)
        lines(z, lower, lty=2, lwd=lwd/2, col=col)
        }
    if (labels[1]==T & length(labels)==1) labels<-seq(along=z)
    if (labels != F) {
        selected<-identify(z, ord.x, labels[good][ord])
        seq(along=x)[good][ord][selected]
        }
    }
    
qq.plot.lm<-function(x, main="", xlab="t Quantiles",
    ylab=paste("Studentized Residuals(",deparse(substitute(x)),")",sep=""),
    line=c("quartiles", "robust"), las=1,
    simulate=F, envelope=.95, labels=F, reps=100, 
    col=palette()[2], lwd=2, pch=1, ...){
    # last modified 1 Feb 2001
    line<-match.arg(line)
    rstudent<-rstudent(x)
    sumry<-summary(x)
    res.df<-sumry$df[2]
    if(!simulate){
        qq.plot.default(rstudent, dist="t", df=res.df-1, line=line,
            main=main, xlab=xlab, ylab=ylab, las=las, envelope=envelope, labels=labels, 
            col=col, lwd=lwd, pch=pch, ...)
        }
    else {
        ord<-order(rstudent)
        ord.x<-rstudent[ord]
        n<-length(ord)
        P<-ppoints(n)
        z<-qt(P, df=res.df-1)
        plot(z, ord.x, xlab=xlab, ylab=ylab, main=main, las=las, pch=pch, col=col)
        yhat<-fitted.values(x)
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
            Q.z<-qt(c(.25,.75),df=res.df-1)
            b<-(Q.x[2]-Q.x[1])/(Q.z[2]-Q.z[1])
            a<-Q.x[1]-b*Q.z[1]
            abline(a, b, col=col, lwd=lwd)
            }
        if (line=="robust"){
            if (!require("MASS", quietly=T)) stop("MASS package not available")
            coef<-coefficients(rlm(ord.x~z))
            a<-coef[1]
            b<-coef[2]
            abline(a, b, col=col, lwd=lwd)
            }
        if (labels[1]==T & length(labels)==1) labels<-seq(along=z)
        if (labels != F) {
            selected<-identify(z, ord.x, labels[ord])
            ord[selected]
            }
        }
    }
 
qq.plot.glm<-function(mod, ...){
    stop("QQ plot for studentized residuals not available for glm")
    }
