# multivariate unconditional Box-Cox transformations (J. Fox)

# last modified 4 June 04 by J. Fox
# (with bug fixes by S. Weisberg)

box.cox.powers<-function(X, start=NULL, hypotheses=NULL, ...){
    modified.power<-function(x, lambda, gm){
        if (lambda == 0) log(x)*gm
        else (gm^(1-lambda))*((x^lambda)-1)/lambda
        }
    neg.kernel.profile.logL<-function(X, lambda, gm){
        for (j in 1:ncol(X)){
            X[,j]<-modified.power(X[,j],lambda[j],gm[j])
            }
        (nrow(X)/2)*log(((nrow(X)-1)/nrow(X))*det(var(X)))
        }
    univ.neg.kernel.logL <- function(x, lambda, gm){
        x <- modified.power(x, lambda, gm)
        (length(x)/2)*log(((length(x)-1)/length(x))*var(x))
        }
    X<-as.matrix(X)
    nc <- ncol(X)
    if(any(X<=0)) stop("All values must be > 0")
    gm<-apply(X, 2, function(x) exp(mean(log(x))))
    if (is.null(start)) {
        start <- rep(1, nc)
        for (j in 1:nc){
            res<- optimize(
                f = function(lambda) univ.neg.kernel.logL(x=X[,j], lambda=lambda, gm=gm[j]),
                lower=-50, upper=+50)
            start[j] <- res$minimum
            }
        }
    res<-optim(start, neg.kernel.profile.logL, hessian=TRUE, method="L-BFGS-B", X=X, gm=gm, ...)
    result<-list()
    result$start<-start
    result$criterion<-res$value
    result$names<-colnames(X)
    result$lambda<-res$par
    result$stderr<-sqrt(diag(inv(res$hessian)))
    result$LR0<-2*(neg.kernel.profile.logL(X,rep(0,nc),gm)-res$value)
    result$LR1<-2*(neg.kernel.profile.logL(X,rep(1,nc),gm)-res$value)
    if (!is.null(hypotheses)) {
        for (i in 1:length(hypotheses)){
            if (length(hypotheses[[i]]) != nc) 
                stop(paste("hypothesis", i, "that powers =", hypotheses[[i]], "does not have", nc, "values"))
            hypotheses[[i]] <- list(test=2*(neg.kernel.profile.logL(X,hypotheses[[i]],gm)-res$value),
                hypothesis=hypotheses[[i]])
            }
        result$hypotheses <- hypotheses
        }
    result$return.code<-res$convergence
    if(result$return.code != 0) 
        warning(paste("Convergence failure: return code =",
            result$return.code))
    class(result)<-"box.cox.powers"
    result
    }
      
print.box.cox.powers<-function(x, digits=4, ...){
    one<-1==length(x$lambda)
    cat(paste("Box-Cox", (if(one) "Transformation to Normality" else "Transformations to Multinormality"),"\n\n"))
    lambda<-x$lambda
    stderr<-x$stderr
    df<-length(lambda)
    result<-cbind(lambda,stderr,lambda/stderr,(lambda-1)/stderr)
    rownames(result)<-x$names
    colnames(result)<-c("Est.Power","Std.Err.",
        "Wald(Power=0)","Wald(Power=1)")
    if (one)rownames(result)<-""
    print(round(result,digits))
    cat(paste("\nL.R. test,", (if(one) "power" else "all powers"), "= 0: ",round(x$LR0,digits),"  df =",df,
        "  p =",round(1-pchisq(x$LR0,df),digits)))
    cat(paste("\nL.R. test,", (if(one) "power" else "all powers"), "= 1: ",round(x$LR1,digits),"  df =",df,
        "  p =",round(1-pchisq(x$LR1,df),digits),"\n"))
    if (!is.null(x$hypotheses)) {
        for (i in 1:length(x$hypotheses)){
            cat(paste("L.R. test, ", (if(one) "power " else "powers "), "= ", 
                paste(x$hypotheses[[i]]$hypothesis,collapse=" "),
                ":  ", round(x$hypotheses[[i]]$test,digits),"   df = ",df,
                "   p = ",round(1-pchisq(x$hypotheses[[i]]$test,df),digits),"\n", sep=""))
            }
        }
    invisible(x)
    }

# summary method retained for backwards compatibility: now identical to print method

summary.box.cox.powers<-function(object, digits=4, ...){
    one<-1==length(object$lambda)
    cat(paste("Box-Cox", (if(one) "Transformation to Normality" else "Transformations to Multinormality"),"\n\n"))
    lambda<-object$lambda
    stderr<-object$stderr
    df<-length(lambda)
    result<-cbind(lambda,stderr,lambda/stderr,(lambda-1)/stderr)
    rownames(result)<-object$names
    colnames(result)<-c("Est.Power","Std.Err.",
        "Wald(Power=0)","Wald(Power=1)")
    if (one)rownames(result)<-""
    print(round(result,digits))
    cat(paste("\nL.R. test,", (if(one) "power" else "all powers"), "= 0: ",round(object$LR0,digits),"  df =",df,
        "  p =",round(1-pchisq(object$LR0,df),digits)))
    cat(paste("\nL.R. test,", (if(one) "power" else "all powers"), "= 1: ",round(object$LR1,digits),"  df =",df,
        "  p =",round(1-pchisq(object$LR1,df),digits),"\n"))
    if (!is.null(object$hypotheses)) {
        for (i in 1:length(object$hypotheses)){
            cat(paste("L.R. test, ", (if(one) "power " else "powers "), "= ", 
                paste(object$hypotheses[[i]]$hypothesis,collapse=" "),
                ":  ", round(object$hypotheses[[i]]$test,digits),"   df = ",df,
                "   p = ",round(1-pchisq(object$hypotheses[[i]]$test,df),digits),"\n", sep=""))
            }
        }
    invisible(object)
    }
