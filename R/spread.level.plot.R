# spread-level plots (J. Fox)

# last modified 2 April 02

slp<-function(x, ...) spread.level.plot(x, ...)

spread.level.plot<-function(x, ...) {
    UseMethod("spread.level.plot")
    }

spread.level.plot.default<-function(x, by, robust.line=any("MASS"==.packages(all=TRUE)), 
        start=0, xlab="Median", ylab="Hinge-Spread", point.labels=TRUE, las=par("las"),
        main=paste("Spread-Level Plot for", deparse(substitute(x)), 
        "by", deparse(substitute(by))), col=palette()[2], pch=1, lwd=2, ...)
    {
    #last modified 23 Feb 2003 by J. Fox
    good<-!(is.na(x) | is.na(by))
    if (sum(good) != length(x)) {
        warning("NAs ignored")
        x<-x[good]
        by<-by[good]
        }    
    min.x<-min(x)
    if (min.x <= -start){
        start<- nice(-min.x +.05*diff(quantile(x,c(.25,.75))), direction="up")
        warning(paste("Start =",start," added to avoid 0 or negative values."))
        }
    if (start !=0) {
        xlab<-paste(xlab, "+", signif(start, 5))
        x<-x+start
        }
    values<-unique(as.character(by))
    result<-matrix(0,length(values),4)
    dimnames(result)<-list(values,c("LowerHinge", "Median", "UpperHinge", "Hinge-Spread"))
    for (i in seq(along=values)){
        five<-fivenum(x[by==values[i]])
        result[i,]<-c(five[2:4],five[4]-five[2])
        }
    medians<-result[,2]
    spreads<-result[,4]
    plot(medians, spreads, log="xy", main=main, xlab=xlab, ylab=ylab, 
        las=las, pch=pch, col=col, ...)
    pos<-ifelse(medians>median(medians), 2, 4)
    if (point.labels) text(medians, spreads, as.character(values), pos=pos, ...)
    if (robust.line){
        if (!require("MASS")) stop("MASS package not available")
        mod<-rlm(log(spreads)~log(medians))
        }
        else mod<-lm(log(spreads)~log(medians), ...)
    ord<-order(medians)
    first<-ord[1]
    last<-ord[length(ord)]
    lines(start+medians[c(first,last)], exp(fitted.values(mod)[c(first,last)]), 
        col=col, lwd=lwd, ...)
    p<-1-(coefficients(mod))[2]
    names(p)<-NULL
    result <- list(Statistics=result[ord,], PowerTransformation=p)
    class(result) <- 'spread.level.plot'
    result
    }
    
spread.level.plot.lm<-function(x, start=0, 
        robust.line=any("MASS"==.packages(all=TRUE)), 
        xlab="Fitted Values",
        ylab="Absolute Studentized Residuals", las=par("las"),
        main=paste("Spread-Level Plot for", deparse(substitute(x))),
        pch=1, col=palette()[2], lwd=2, ...)
    {
    #last modified 23 Feb 2003 by J. Fox
    resid<-na.omit(abs(rstudent(x)))
    fitval<-na.omit(fitted.values(x))
    min<-min(fitval)
    if (min <= -start) {
        start<- nice(-min +.05*diff(quantile(fitval,c(.25,.75))), direction='up')
        warning(paste("Start = ", start, 
            "added to fitted values to avoid 0 or negative values."))
        }
    if (start !=0) xlab<-paste(xlab, "+", signif(start, 5))
    plot(fitval+start, resid, log="xy", main=main, xlab=xlab, ylab=ylab, 
        las=las, col=col, pch=pch, ...)
    if (robust.line){
        if (!require("MASS")) stop("MASS package not available")
        mod<-rlm(log(resid)~log(fitval+start))
        }
        else mod<-lm(log(resid)~log(fitval+start), ...)
    first<-which.min(fitval) 
    last<-which.max(fitval) 
    lines((fitval+start)[c(first,last)], exp(fitted.values(mod)[c(first,last)]), 
        lwd=lwd, col=col, ...)
    p<-1-(coefficients(mod))[2]
    names(p)<-NULL
    result <- list(PowerTransformation=p)
    class(result) <- 'spread.level.plot'
    result
    }
  
spread.level.plot.formula<-function (formula, data=NULL, subset, na.action, 
    main=paste("Spread-Level Plot for", varnames[response], "by", varnames[-response]), ...) {
    if (missing(na.action)) 
        na.action <- options()$na.action
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
        m$data <- as.data.frame(data)
    m$... <- m$main <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, sys.frame(sys.parent()))
    response <- attr(attr(mf, "terms"), "response")
    varnames <- names(mf)
    if (!response) stop ("No response variable specified")
    if (length(varnames)>2) stop("Right-hand side of model has more than one variable")
    x <- mf[[response]]
    by <- mf[[varnames[-response]]]
    spread.level.plot.default(x, by, main=main, ...)
    }
    
print.spread.level.plot <- function(x, ...){
    if (!is.null(x$Statistics)) print(x$Statistics, ...)
    cat('\nSuggested power transformation: ', x$PowerTransformation,'\n')
    invisible(x)
    }
