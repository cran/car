# Linear hypothesis tests for lm and glm (J. Fox)

lht<-function(...) linear.hypothesis(...)

linear.hypothesis<-function (model, ...) {
    UseMethod("linear.hypothesis")
    }
    
linear.hypothesis.lm<-function(model, hypothesis.matrix, rhs=0, 
        summary.model=summary(model, corr = FALSE),
        white.adjust=F, error.SS, error.df) {
    # last modified by J.Fox 13 Dec 2000
    if (is.aliased(model)) stop ("One or more terms aliased in model.")
    s2<-if (missing(error.SS)) summary.model$sigma^2
        else error.SS/error.df
    V<-if (white.adjust==FALSE) summary.model$cov.unscaled
        else hccm(model, type=white.adjust)/s2
    b<-coefficients(model)
    L<-if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
        else hypothesis.matrix
    q<-nrow(L)
    SSH<-t(L %*% b - rhs) %*% inv(L %*% V %*% t(L)) %*% (L %*% b - rhs)
    F<-SSH/(q*s2)
    df<-if (missing(error.df)) model$df.residual
        else error.df
    p<-1-pf(F, q, df)
    result<-list(SSH=SSH[1,1], SSE=s2*df, F=F[1,1], Df=c(q, df), p=p[1,1])
    class(result)<-"F.test"
    result
    }

linear.hypothesis.glm<-function(model, hypothesis.matrix, rhs=0, 
        summary.model=summary(model, corr = FALSE)) {
    # last modified by J.Fox 13 Dec 2000
    if (is.aliased(model)) stop ("One or more terms aliased in model.")
    V<-summary.model$dispersion*summary.model$cov.unscaled
    b<-coefficients(model)
    L<-if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
        else hypothesis.matrix
    q<-nrow(L)
    Wald<-t(L %*% b - rhs) %*% inv(L %*% V %*% t(L)) %*% (L %*% b - rhs)
    p<-1-pchisq(Wald, q)
    result<-list(test="Wald Test",ChiSquare=Wald[1,1], Df=q, p=p[1,1])
    class(result)<-"chisq.test"
    result
    }
 

 print.chisq.test<-function(x){
    title<-if (!is.null(x$test)) x$test else "Chisquare Test"
    cat(title,"\n")
    if (!is.null(x$formula)) cat(x$formula.name, 
        "formula:", as.character(x$formula), "\n")
    cat("Chisquare =", x$ChiSquare,"   Df =", x$Df,
        "    p =", x$p, "\n")
    invisible(x)
    }
 
 print.F.test<-function(x){
    title<-if (!is.null(x$test)) x$test else "F-Test"
    cat(title,"\n")
    if (!is.null(x$formula)) cat(x$formula.name, 
        "formula:", as.character(x$formula), "\n")
    cat("SS =", x$SSH, "    SSE =", x$SSE, "    F =", x$F,
        " Df =", x$Df[1], "and", x$Df[2], "    p =", x$p, "\n")
    invisible(x)
    }
