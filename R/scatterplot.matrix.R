# fancy scatterplot matrices (J. Fox)

# last modified: 24 August 2009 by J. Fox

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
        groups <- FALSE
        }
    else{
        groups <- TRUE
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

scatterplot.matrix.default<-function(x, labels=colnames(x), 
    diagonal=c("density", "boxplot", "histogram", "oned", "qqplot", "none"), adjust=1, nclass,
    plot.points=TRUE, smooth=TRUE, span=.5, reg.line=lm, transform=FALSE,
    ellipse=FALSE, levels=c(.5, .9), robust=FALSE,
    groups=FALSE, by.groups=FALSE,
    col=palette(), pch=1:n.groups, lwd=1, lwd.smooth=lwd,
    cex=par("cex"), cex.axis=par("cex.axis"), cex.labels=NULL, 
    cex.main=par("cex.main"),
    legend.plot=length(levels(groups)) > 1, ...){
    # last modified 16 Jan 2005 by J. Fox
    if (groups[1] != FALSE){
        x<-na.omit(cbind(as.data.frame(groups),x))
        groups<-as.factor(as.character(x[,1]))
        x<-x[,-1]
        }
        else x<-na.omit(x)
    if (missing(nclass)) nclass<-n.bins(x[,1])
    reg<-function(x, y, col){
        mod<-reg.line(y~x)
        y.hat<-fitted.values(mod)
        x<-model.matrix(mod)[,2]
        min<-which.min(x)
        max<-which.max(x)
        lines(c(x[min],x[max]),c(y.hat[min],y.hat[max]), lty=2, lwd=lwd, col=col)
        }
    # The following panel function adapted from Richard Heiberger
    panel.density<-function(x, ...){
        dens.x <- density(x, adjust = adjust)
        lines(dens.x$x, min(x) + dens.x$y * diff(range(x))/diff(range(dens.x$y)))
        points(x, rep(min(x), length(x)), pch = "|", col = col[1])
        }
    panel.histogram<-function(x, ...){
        par(new=TRUE)
        hist(x, main="", axes=FALSE, nclass=nclass, col=col[2])
        }
    panel.boxplot<-function(x, ...){
        par(new=TRUE)
        boxplot(x, axes=FALSE, main="", col=col[2])
        }
    # The following panel function adapted from Richard Heiberger
    panel.oned <- function(x, ...) {
      range <- range(x)
      delta <- diff(range)/50
      y <- mean(range)
      segments(x-delta, x, x+delta, x, col = col[1])
    }
    panel.qqplot<-function(x, ...){
        par(new=TRUE)
        qqnorm(x, axes=FALSE, xlab="", ylab="", main="", col=col[1])
        qqline(x)
        }
    panel.blank<-function(x, ...) NULL
    which.fn<-match(match.arg(diagonal),
        c("density", "boxplot", "histogram", "oned", "qqplot", "none"))
    diag<-list(panel.density, panel.boxplot, panel.histogram, panel.oned, panel.qqplot, panel.blank)[[which.fn]]
    groups<-as.factor(if(FALSE==groups[1]) rep(1, length(x[,1])) else groups)
    n.groups<-length(levels(groups))
    if (n.groups >= length(col)) stop("number of groups exceeds number of available colors")
    if (transform != FALSE | length(transform) == ncol(x)){
        if (transform == TRUE & length(transform) == 1) transform <- box.cox.powers(x)$lambda
        for (i in 1:ncol(x)){
            x[,i]<-box.cox(x[,i], transform[i])
            labels[i] <- paste(labels[i], "^(", round(transform[i],2), ")", sep="")
            }
        }          
    pairs(x, labels=labels,
        cex.axis=cex.axis, cex.main=cex.main, cex.labels=cex.labels, cex=cex,
        diag.panel=diag,
        panel=function(x, y, ...){ 
            for (i in 1:n.groups){
                subs<-groups==levels(groups)[i]
                if (plot.points) points(x[subs], y[subs], pch=pch[i], col=col[i+1], cex=cex)
                if (smooth & by.groups) lines(lowess(x[subs], y[subs]), col=col[i+1], lwd=lwd.smooth)
                if (is.function(reg.line) & by.groups) reg(x[subs], y[subs], col=col[i+1])
                if (ellipse  & by.groups) data.ellipse(x[subs], y[subs], plot.points=FALSE, 
                    levels=levels, col=col[i+1], robust=robust, lwd=1)
                }
            if (!by.groups){
                if (is.function(reg.line)) abline(reg.line(y~x),lty=2, lwd=lwd, col=col[1])
                if (smooth) lines(lowess(x,y, f=span), lwd=lwd.smooth, col=col[1])
                if (ellipse) data.ellipse(x, y, plot.points=FALSE, levels=levels, col=col[1],
                    robust=robust, lwd=1)
                }
            }, ...
        )
    if(legend.plot) {
        frac<-1/ncol(x)
        legend(1 - .95*frac, 0.8*frac,
            legend=levels(groups), pch=pch, col=col[2:(n.groups+1)],
            cex=cumprod(par("fin"))[2]*sqrt(frac)/(sqrt(n.groups)*20))
        }
    }

spm<-function(x, ...){
    scatterplot.matrix(x, ...)
    }            
