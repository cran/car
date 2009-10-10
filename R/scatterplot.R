# fancy scatterplots  (J. Fox)

# last modified 24 August 2009

scatterplot<-function(x, ...){
    # last modified 28 Jan 2001 by J. Fox
    UseMethod("scatterplot", x)
    }
    
scatterplot.formula<-function (formula, data, xlab, ylab, legend.title, subset, labels=FALSE, ...) {
    # last modified 6 Jan 2004 by J. Fox
    na.save <- options(na.action=na.omit)
    on.exit(options(na.save))
    na.pass<-function(dframe) dframe
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
        m$data <- as.data.frame(data)
    m$na.action <- na.pass
    m$labels <- NULL
    m$... <- NULL
    m$xlab <- m$ylab <- NULL
    m[[1]] <- as.name("model.frame")
    if (!inherits(formula, "formula") | length(formula) != 3) 
        stop("invalid formula")    
    formula<-as.character(c(formula))
    formula<-as.formula(sub("\\|", "+", formula))
    m$formula<-formula
    if (missing(data)){ 
        X <- na.omit(eval(m, parent.frame()))
        if (labels[1] != FALSE) labels<-labels[as.numeric(gsub("X","", row.names(X)))]
        }
    else{
        if (labels[1] != FALSE) row.names(data)<-labels
        X <- eval(m, parent.frame())
        if (labels[1] != FALSE) labels<-row.names(X)
        }
    names<-names(X)
    if (missing(xlab)) xlab<-names[2]
    if (missing(ylab)) ylab<-names[1]
    if (ncol(X) == 2) scatterplot(X[,2], X[,1],  xlab=xlab, ylab=ylab, 
            labels=labels, ...)
    else {
        if (missing(legend.title)) legend.title <- names[3]
        scatterplot(X[,2], X[,1], groups=X[,3], xlab=xlab, ylab=ylab,  
                legend.title=legend.title, labels=labels, ...)
        }
    }


scatterplot.default<-function(x, y, smooth=TRUE, span=.5, reg.line=lm, boxplots="xy",
    xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), las=par("las"),
    lwd=1, lwd.smooth=lwd, labels=FALSE, log="", jitter=list(), xlim=NULL, ylim=NULL,
    cex=par("cex"), cex.axis=par("cex.axis"), cex.lab=par("cex.lab"), 
    cex.main=par("cex.main"), cex.sub=par("cex.sub"),
    groups=FALSE, by.groups=!(groups[1]==FALSE), legend.title=deparse(substitute(groups)), 
    ellipse=FALSE, levels=c(.5, .9), robust=FALSE,
    col=rep(palette(), length.out=n.groups + 1), pch=1:n.groups, 
    legend.plot=length(levels(groups)) > 1, reset.par=TRUE, ...){
    # last modified 21 Oct 2008 by J. Fox
    lowess.line<-function(x, y, col) {
        x<-if (0==length(grep("x", log))) x else log(x)
        y<-if (0==length(grep("y", log))) y else log(y)
        valid<-!(is.na(x) | is.na(y))
        fit<-lowess(x[valid],y[valid],f=span)
        x<-if (0==length(grep("x", log))) fit$x else exp(fit$x)
        y<-if (0==length(grep("y", log))) fit$y else exp(fit$y)
        lines(x, y, lwd=lwd.smooth, col=col)
        }
    reg<-function(x, y, col){
        x<-if (0==length(grep("x", log))) x else log(x)
        y<-if (0==length(grep("y", log))) y else log(y)
        mod<-reg.line(y~x)
        y.hat<-fitted.values(mod)
        x<-model.matrix(mod)[,2]
        min<-which.min(x)
        max<-which.max(x)
        if (0==length(grep("x", log))){
            x1<-x[min]
            x2<-x[max]
            }
        else {
            x1<-exp(x[min])
            x2<-exp(x[max])
            }
        if (0==length(grep("y", log))){
            y1<-y.hat[min]
            y2<-y.hat[max]
            }
        else {
            y1<-exp(y.hat[min])
            y2<-exp(y.hat[max])
            }
        lines(c(x1,x2),c(y1,y2), lty=2, lwd=lwd, col=col)
        }
    hbox<-function(x){
        if (length(grep("x", log))==0){
            log.x<-""
            .x<-x
            }
        else {
            log.x<-"x"
            .x<-log(x)
            }
        plot(x, seq(0,1,length=length(x)), type="n", axes=FALSE, xlab="", ylab="", log=log.x, xlim=xlim)
        res<-boxplot.stats(.x, coef = 1.5, do.conf=FALSE)
        if (length(grep("x", log))!=0){
            res$stats<-exp(res$stats)
            if (!is.null(res$out)) res$out<-exp(res$out)
            }
        LW<-res$stats[1]
        Q1<-res$stats[2]
        M<-res$stats[3]
        Q3<-res$stats[4]
        UW<-res$stats[5]
        lines(c(Q1,Q1,Q3,Q3,Q1),c(0,1,1,0,0))
        lines(c(M,M),c(0,1))
        lines(c(LW,Q1),c(.5,.5))
        lines(c(Q3,UW),c(.5,.5))
        if (!is.null(res$out)) points(res$out,rep(.5, length(res$out)), cex=cex)
        }
    vbox<-function(y){
        if (length(grep("y", log))==0){
            log.y<-""
            .y<-y
            }
        else {
            log.y<-"y"
            .y<-log(y)
            }
        plot(seq(0,1,length=length(y)), y, type="n", axes=FALSE, xlab="", ylab="", log=log.y, ylim=ylim)
        res<-boxplot.stats(.y, coef = 1.5, do.conf=FALSE)
        if (length(grep("y", log))!=0){
            res$stats<-exp(res$stats)
            if (!is.null(res$out)) res$out<-exp(res$out)
            }
        LW<-res$stats[1]
        Q1<-res$stats[2]
        M<-res$stats[3]
        Q3<-res$stats[4]
        UW<-res$stats[5]
        lines(c(0,1,1,0,0),c(Q1,Q1,Q3,Q3,Q1))
        lines(c(0,1),c(M,M))
        lines(c(.5,.5),c(LW,Q1))
        lines(c(.5,.5),c(Q3,UW))
        if (!is.null(res$out)) points(rep(.5, length(res$out)),res$out, cex=cex)
        }
    mar<-par("mar")
    mfcol<-par("mfcol")
    if (reset.par) on.exit(par(mar=mar, mfcol=mfcol))
    if(FALSE==boxplots) boxplots<-""
    if (groups[1] != FALSE){
        legend.title # force evaluation
        if (labels[1] != FALSE){
            data<-na.omit(data.frame(groups,x,y,labels))
            groups<-data[,1]
            .x<-data[,2]
            .y<-data[,3]
            labels<-data[,4]
            }
        else {
            data<-na.omit(data.frame(groups,x,y))
            groups<-data[,1]
            .x<-data[,2]
            .y<-data[,3]
            }
			top <- if (legend.plot) 4 + length(levels(as.factor(groups)))
				else mar[3]
        }
    else{
        .x<-x
        .y<-y
        top <- mar[3]
        }
    groups<-as.factor(if(FALSE == groups[1]) rep(1, length(.x)) else as.character(groups))
    layout(matrix(c(1,0,3,2),2,2),
        widths = c(5,95),
        heights= c(95,5))
    par(mar=c(mar[1],0,top,0))
    if (length(grep("y",boxplots))>0) vbox(.y) else plot(0,0,xlab="",ylab="",axes=FALSE,type="n", xlim=xlim, ylim=ylim)
    par(mar=c(0,mar[2],0,mar[4]))
    if (length(grep("x",boxplots))>0) hbox(.x) else plot(0,0,xlab="",ylab="",axes=FALSE,type="n", xlim=xlim, ylim=ylim)
    par(mar=c(mar[1:2], top, mar[4]))
    plot(.x, .y, xlab=xlab, ylab=ylab, las=las, log=log, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab,
        cex.main=cex.main, cex.sub=cex.sub, type="n", xlim=xlim, ylim=ylim, ...)
    n.groups<-length(levels(groups))
    if (n.groups >= length(col)) stop("number of groups exceeds number of available colors")
    for (i in 1:n.groups){
        subs<-groups==levels(groups)[i]
        points(if (is.null(jitter$x) || jitter$x == 0) .x[subs] else jitter(.x[subs], factor=jitter$x), 
            if (is.null(jitter$y) || jitter$y == 0) .y[subs] else jitter(.y[subs], factor=jitter$y), 
            pch=pch[i], col=col[i+1], cex=cex)
        if (smooth & by.groups) lowess.line(.x[subs], .y[subs], col=col[i+1])
        if (is.function(reg.line) & by.groups) reg(.x[subs], .y[subs], col=col[i+1])
        if (ellipse  & by.groups) data.ellipse(.x[subs], .y[subs], plot.points=FALSE, 
            levels=levels, col=col[i+1], robust=robust)
        }
    if (!by.groups){
        if (smooth) lowess.line(.x, .y, col=col[1])
        if (is.function(reg.line)) reg(.x, .y, col=col[1])
        if (ellipse) data.ellipse(.x, .y, plot.points=FALSE, levels=levels, col=col[1],
            robust=robust)
        }
    if(legend.plot) {
        xpd <- par(xpd=TRUE)
        on.exit(par(xpd=xpd), add=TRUE)
        usr <- par("usr")
        legend(usr[1], usr[4] + 1.2*top*strheight("x"), legend=levels(groups), 
            pch=pch, col=col[2:(n.groups+1)], pt.cex=cex, cex=cex.lab, title=legend.title)
        }
    if (labels[1]==TRUE & length(labels)==1) labels<-seq(along=.x)
    indices<-if (labels[1] != FALSE) identify(.x, .y, labels)
    if (is.null(indices)) invisible(indices) else indices
    }

sp<-function(...) scatterplot(...)
