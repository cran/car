# Box-Tidwell transformations (J. Fox)

# last modified 2 April 02 by J. Fox

box.tidwell<-function(y, ...){
    UseMethod("box.tidwell")
    }

box.tidwell.formula<-function(formula, other.x=NULL, data=NULL, subset, na.action=options()$na.action, 
    verbose=FALSE, tol=.001, max.iter=25, ...) {
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
        m$data <- as.data.frame(data)
    m$formula<-if (is.null(other.x)) formula
        else as.formula(paste(formula[2], "~", formula[3], "+", other.x[2]))
    m$max.iter<-m$tol<-m$verbose<-m$family<-m$other.x <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, sys.frame(sys.parent()))
    response <- attr(attr(mf, "terms"), "response")
    if (!response) stop(paste("No response variable in model."))
    X1<-model.matrix(formula, data=mf)[,-1]
    X2<-if (is.null(other.x)) NULL
        else model.matrix(other.x, data=mf)[,-1]
    y<-model.response(mf, "numeric")
    box.tidwell.default(y, X1, X2, max.iter=max.iter, tol=tol, verbose=verbose, ...)
    }

box.tidwell.default<-function(y, x1, x2=NULL, max.iter=25, tol=.001, verbose=FALSE, ...) {
    # last modified 19 Sept 02 by J. Fox
    x1<-as.matrix(x1)
    var.names<-if(is.null(colnames(x1))) 1:ncol(x1) else colnames(x1)
    k.x1<-length(var.names)
    x.log.x<-x1*log(x1)
    mod.1<-lm(y~cbind(x1, x2), ...)
    mod.2<-lm(y~cbind(x.log.x, x1, x2), ...)
    sumry<-summary(mod.2)
    seb<-sqrt(diag(vcov(mod.2)))
    t.vals<-((coefficients(mod.2))/seb)[2:(1+k.x1)]
    initial<-powers<-1+coefficients(mod.2)[2:(1+k.x1)]/coefficients(mod.1)[2:(1+k.x1)]
    pvalues<-2*(1-pnorm(abs(t.vals)))
    iter<-0
    last.powers<-1
    while ( (max(abs((powers-last.powers)/(powers+tol))) > tol) &
        (iter <= max.iter) ) {
        iter<-iter+1
        x1.p<-x1^matrix(powers, nrow=nrow(x1), ncol=ncol(x1), byrow=TRUE)
        x.log.x<-x1.p*log(x1.p)
        mod.1<-lm(y~cbind(x1.p, x2), ...)
        mod.2<-lm(y~cbind(x.log.x, x1.p, x2), ...)
        last.powers<-powers
        powers<-powers * 
            (1+coefficients(mod.2)[2:(1+k.x1)]/coefficients(mod.1)[2:(1+k.x1)])
        if (verbose) cat(" iter =", iter, "    powers =", powers, "\n")
        }
    if (iter > max.iter) warning("maximum iterations exceeded")
    result<-rbind(initial,t.vals, pvalues, powers)
    rownames(result)<-c("Initial Power","Score Statistic","p-value","MLE of Power")
    colnames(result)<-names(powers)<-var.names
    result<-list(result=result, iterations=iter)
    class(result)<-"box.tidwell"
    result
    }

print.box.tidwell<-function(x, digits=5, ...){ 
    # last modified 15 Dec 2000 by J. Fox  
    print(round(x$result, digits))
    cat("\niterations = ", x$iterations,"\n")
    }
  
