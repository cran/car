# fancy scatterplot matrices (J. Fox)

scatterplot.matrix<-function(x, ...){
    UseMethod("scatterplot.matrix")
    }

scatterplot.matrix.formula<-function (formula, data=NULL, subset,  ...) {
    # last modified 1 Feb 2001 by J. Fox
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m$formula<-NULL
    if (!inherits(formula, "formula") | length(formula) != 2) 
        stop("invalid formula")
    rhs <- formula[[2]]
    if ("|"!=deparse(rhs[[1]])){
        groups <- F
        }
    else{
        groups <- T
        formula<-as.character(c(formula))
        formula<-as.formula(sub("\\|", "+", formula))   
        }
    m$formula<-formula
    X <- eval(m, sys.frame(sys.parent()))
    if (!groups) scatterplot.matrix(X, ...)
        else{
        ncol<-ncol(X)
        scatterplot.matrix.default(X[,-ncol], groups=X[,ncol], ...)
        }
    }

scatterplot.matrix.default<-function(data, labels=colnames(data), 
    diagonal=c("density", "boxplot", "histogram", "qqplot", "none"), adjust=1, nclass,
    plot.points=TRUE, smooth=TRUE, span=.5, reg.line=lm, transform=FALSE,
    ellipse=FALSE, levels=c(.5, .9), robust=FALSE,
    groups=FALSE, by.groups=FALSE,
    col=palette(), pch=1:n.groups, lwd=1,
    legend.plot=length(levels(groups)) > 1){
    # last modified 23 April 2001 by J. Fox
    if (groups != FALSE){
        data<-na.omit(cbind(as.data.frame(groups),data))
        groups<-as.factor(as.character(data[,1]))
        data<-data[,-1]
        }
    if (missing(nclass)) nclass<-n.bins(data[,1])
    reg<-function(x, y, col){
        mod<-reg.line(y~x)
        y.hat<-fitted.values(mod)
        x<-model.matrix(mod)[,2]
        min<-which.min(x)
        max<-which.max(x)
        lines(c(x[min],x[max]),c(y.hat[min],y.hat[max]), lty=2, lwd=lwd, col=col)
        }
    panel.density<-function(x){
        par(new=T)
        plot(density(x, adjust=adjust), axes=F, main="")
        points(x, rep(0,length(x)), pch="|", col=col[1])
        }
    panel.histogram<-function(x){
        par(new=T)
        hist(x, main="", axes=F, nclass=nclass, col=col[1])
        }
    panel.boxplot<-function(x){
        par(new=T)
        boxplot(x, axes=F, main="", col=col[1])
        }
    panel.qqplot<-function(x){
        par(new=T)
        qqnorm(x, axes=F, xlab="", ylab="", main="", col=col[1])
        qqline(x)
        }
    panel.blank<-function(x) NULL
    which.fn<-match(match.arg(diagonal), c("density", "boxplot", "histogram", "qqplot", "none"))
    diag<-list(panel.density, panel.boxplot, panel.histogram, panel.qqplot, panel.blank)[[which.fn]]
    groups<-as.factor(if(FALSE==groups) rep(1, length(data[,1])) else groups)
    n.groups<-length(levels(groups))
    if (n.groups >= length(col)) stop("number of groups exceeds number of available colors")
    if (transform != F | length(transform) == ncol(data)){
        if (transform == T & length(transform) == 1) transform <- box.cox.powers(data)$lambda
        for (i in 1:ncol(data)){
            data[,i]<-box.cox(data[,i], transform[i])
            labels[i] <- paste(labels[i], "^(", round(transform[i],2), ")", sep="")
            }
        }          
    pairs(data, labels=labels,
        diag.panel=diag,
        panel=function(x, y, ...){ 
            for (i in 1:n.groups){
                subs<-groups==levels(groups)[i]
                if (plot.points) points(x[subs], y[subs], pch=pch[i], col=col[i+1])
                if (smooth & by.groups) lines(lowess(x[subs], y[subs]), col=col[i+1])
                if (is.function(reg.line) & by.groups) reg(x[subs], y[subs], col=col[i+1])
                if (ellipse  & by.groups) data.ellipse(x[subs], y[subs], plot.points=F, 
                    levels=levels, col=col[i+1], robust=robust)
                }
            if (!by.groups){
                if (is.function(reg.line)) abline(reg.line(y~x),lty=2, lwd=lwd, col=col[1])
                if (smooth) lines(lowess(x,y, f=span), lwd=lwd, col=col[1])
                if (ellipse) data.ellipse(x, y, plot.points=F, levels=levels, col=col[1],
                    robust=robust)
                }
            }
        )
    if(legend.plot) {
        frac<-1/ncol(data)
        legend(1 - .95*frac, 0.8*frac,
            legend=levels(groups), pch=pch, col=col[2:(n.groups+1)], 
            cex=cumprod(par("fin"))[2]*sqrt(frac)/(sqrt(n.groups)*20))
        }
    }

spm<-function(x, ...){
    scatterplot.matrix(x, ...)
    }            
