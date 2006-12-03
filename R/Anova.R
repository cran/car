# Type II and III tests for linear, generalized linear, and other models (J. Fox)

# last modified 2 December 06

relatives<-function(term, names, factors){
    is.relative<-function(term1, term2) {
        all(!(factors[,term1]&(!factors[,term2])))
        }
    if(length(names)==1) return(NULL)
    which.term<-which(term==names)
    (1:length(names))[-which.term][sapply(names[-which.term], function(term2) is.relative(term, term2))]
    }
    

Anova<-function(mod, ...){
    UseMethod("Anova", mod)
    }

 # linear models
 
Anova.lm<-function(mod, error, type=c("II","III", 2, 3), ...){
    type<-match.arg(type)
    switch(type,
        II=Anova.II.lm(mod, error, ...),
        III=Anova.III.lm(mod, error, ...),
        "2"=Anova.II.lm(mod, error, ...),
        "3"=Anova.III.lm(mod, error, ...))
    }

Anova.aov <- function(mod, ...){
    # last modified 8 Mar 2002 by J. Fox
    class(mod) <- "lm"
    Anova.lm(mod, ...)
    }

        # type II
        
Anova.II.lm<-function(mod, error, ...){
    # last modified by J.Fox 04 April 2005
    if (!missing(error)){
        sumry<-summary(error, corr=FALSE)
        s2<-sumry$sigma^2
        error.df<-error$df.residual
        error.SS<-s2*error.df
        }
    SS.term<-function(term){
        which.term<-which(term==names)
        subs.term<-which(assign==which.term)
        relatives<-relatives(term, names, fac)
        subs.relatives<-NULL
        for (relative in relatives) subs.relatives<-c(subs.relatives, which(assign==relative))
        hyp.matrix.1<-I.p[subs.relatives,]
        hyp.matrix.2<-I.p[c(subs.relatives,subs.term),]
        SS1<-if (length(subs.relatives)==0) 0 
            else linear.hypothesis(mod, hyp.matrix.1, summary.model=sumry, ...)$"Sum of Sq"[2]
        SS2<-linear.hypothesis(mod, hyp.matrix.2, summary.model=sumry, ...)$"Sum of Sq"[2]
        SS1 - SS2
        }
    fac<-attr(mod$terms, "factors")
    intercept<-has.intercept(mod)
    p<-length(coefficients(mod))
    I.p<-diag(p)
    assign<-mod$assign
    names<-term.names(mod)
    if (intercept) names<-names[-1]
    n.terms<-length(names)
    SS<-rep(0, n.terms+1)
    df<-rep(0, n.terms+1)
    f<-rep(0, n.terms+1)
    p<-rep(0, n.terms+1)
    sumry<-summary(mod, corr = FALSE)
    SS[n.terms+1]<-if (missing(error)) sumry$sigma^2*mod$df.residual else error.SS   
    df[n.terms+1]<-if (missing(error)) mod$df.residual else error.df
    f[n.terms+1]<-NA
    p[n.terms+1]<-NA
    for (i in 1:n.terms){
        SS[i]<-SS.term(names[i])
        df[i]<-df.terms(mod, names[i])
        f[i]<-df[n.terms+1]*SS[i]/(df[i]*SS[n.terms+1])
        p[i]<-1-pf(f[i],df[i],df[n.terms+1])
        }    
    result<-data.frame(SS, df, f, p)
    row.names(result)<-c(names,"Residuals")
    names(result)<-c("Sum Sq", "Df", "F value", "Pr(>F)")
    class(result)<-c("anova","data.frame")
    attr(result,"heading")<-c("Anova Table (Type II tests)\n", paste("Response:", responseName(mod)))
    result
    }

        # type III
        
Anova.III.lm<-function(mod, error, ...){
    # last modified by J.Fox 4 April 2005
    if (!missing(error)){
        sumry<-summary(error, corr=FALSE)
        s2<-sumry$sigma^2
        error.df<-error$df.residual
        error.SS<-s2*error.df
        }
    intercept<-has.intercept(mod)
    p<-length(coefficients(mod))
    I.p<-diag(p)
    Source<-term.names(mod)
    n.terms<-length(Source)
    SS<-rep(0, n.terms+1)
    df<-rep(0, n.terms+1)
    f<-rep(0, n.terms+1)
    p<-rep(0, n.terms+1)
    assign<-mod$assign
    sumry<-summary(mod, corr = FALSE)
    for (term in 1:n.terms){
        subs<-which(assign==term-intercept)
        hyp.matrix<-I.p[subs,]
        test<-if (missing(error)) linear.hypothesis(mod, hyp.matrix, summary.model=sumry, ...)
            else linear.hypothesis(mod, hyp.matrix, error.SS=error.SS, error.df=error.df, 
                summary.model=sumry, ...)
        SS[term]<- -test$"Sum of Sq"[2]
        df[term]<- -test$"Df"[2]
        f[term]<-test$"F"[2]
        p[term]<-test$"Pr(>F)"[2]
        }
     Source[n.terms+1]<-"Residuals"
     df.res<-if (missing(error)) mod$df.residual
        else error.df     
     s2<-sumry$sigma^2
     SS[n.terms+1]<-if (missing(error)) s2*df.res
        else error.SS
     df[n.terms+1]<-df.res
     f[n.terms+1]<-NA
     p[n.terms+1]<-NA
     result<-data.frame(SS, df, f, p)
     row.names(result)<-Source
     names(result)<-c("Sum Sq", "Df", "F value", "Pr(>F)")
     class(result)<-c("anova","data.frame")
     attr(result,"heading")<-c("Anova Table (Type III tests)\n", paste("Response:", responseName(mod)))
     result
     }

    # generalized linear models
    
Anova.glm<-function(mod, type=c("II","III", 2, 3), test.statistic=c("LR", "Wald", "F"), 
    error, error.estimate=c("pearson", "dispersion", "deviance"), ...){
    #last modified by J. Fox 1 Dec 2006
    type<-match.arg(type)
    test.statistic<-match.arg(test.statistic)
    error.estimate<-match.arg(error.estimate)
    switch(type,
        II=switch(test.statistic,
            LR=Anova.II.LR.glm(mod),
            Wald=Anova.II.Wald.glm(mod),
            F=Anova.II.F.glm(mod, error, error.estimate)),
        III=switch(test.statistic,
            LR=Anova.III.LR.glm(mod),
            Wald=Anova.III.Wald.glm(mod),
            F=Anova.III.F.glm(mod, error, error.estimate)),
        "2"=switch(test.statistic,
            LR=Anova.II.LR.glm(mod),
            Wald=Anova.II.Wald.glm(mod),
            F=Anova.II.F.glm(mod, error, error.estimate)),
        "3"=switch(test.statistic,
            LR=Anova.III.LR.glm(mod),
            Wald=Anova.III.Wald.glm(mod),
            F=Anova.III.F.glm(mod, error, error.estimate)))
    }

    
        # type III
        
            # Wald test
        
Anova.III.Wald.glm<-function(mod, ...){
    # last modified by J.Fox 11 Dec 2000
    intercept<-has.intercept(mod)
    p<-length(coefficients(mod))
    I.p<-diag(p)
    Source<-term.names(mod)
    n.terms<-length(Source)
    Wald<-rep(0, n.terms)
    df<-rep(0, n.terms)
    p<-rep(0, n.terms)
    assign<-attr(model.matrix(mod),"assign")
    sumry<-summary(mod, corr=FALSE)
    for (term in 1:n.terms){
        subs<-which(assign==term-intercept)
        hyp.matrix<-I.p[subs,]
        test<-linear.hypothesis(mod, hyp.matrix, summary.model=sumry)
        Wald[term]<-test$Chisq[2]
        df[term]<- -test$Df[2]
        p[term]<- test$"Pr(>Chisq)"[2]
        }
     result<-data.frame(Wald, df, p)
     row.names(result)<-Source
     names(result)<-c("Wald Chisq","Df","Pr(>Chisq)")
     class(result)<-c("anova","data.frame")
     attr(result,"heading")<-c("Anova Table (Type III tests)\n", paste("Response:", responseName(mod)))
     result
     }
     
            # LR test

Anova.III.LR.glm<-function(mod, ...){
    Source<-if (has.intercept(mod)) term.names(mod)[-1]
        else term.names(mod)
    n.terms<-length(Source)
    LR<-rep(0, n.terms)
    df<-rep(0, n.terms)
    p<-rep(0, n.terms)
    dispersion<-summary(mod, corr = FALSE)$dispersion
    deviance<-deviance(mod)/dispersion
    for (term in 1:n.terms){
        mod.1<-drop1(mod, scope=
            eval(parse(text=paste("~",Source[term]))))
        LR[term]<-(mod.1$Deviance[2]/dispersion)-deviance
        df[term]<-mod.1$Df[2]
        p[term]<-1-pchisq(LR[term], df[term])
        }
     result<-data.frame(LR, df, p)
     row.names(result)<-Source
     names(result)<-c("LR Chisq","Df","Pr(>Chisq)")
     class(result)<-c("anova","data.frame")
     attr(result,"heading")<-c("Anova Table (Type III tests)\n", paste("Response:", responseName(mod)))
     result
     }

            # F test

Anova.III.F.glm<-function(mod, error, error.estimate, ...){
    # last modified by J. Fox 25 Apr 2003
    fam <- family(mod)$family
    if (fam == "binomial" || fam == "poisson") 
        warning("dispersion parameter estimated from the Pearson residuals, not taken as 1")
    if (missing(error)) error<-mod
    df.res <- df.residual(error)
    error.SS<-switch(error.estimate,
        pearson=sum(residuals(error, "pearson")^2),
        dispersion=df.res*summary(error, corr = FALSE)$dispersion,
        deviance=deviance(error))
    Source<-if (has.intercept(mod)) term.names(mod)[-1]
        else term.names(mod)
    n.terms<-length(Source)
    p <- df <- f <- SS <-rep(0, n.terms+1)
    f[n.terms+1] <- p[n.terms+1] <- NA
    df[n.terms+1]<-df.res
    SS[n.terms+1]<-error.SS
    dispersion<-error.SS/df.res
    deviance<-deviance(mod)
    for (term in 1:n.terms){
        mod.1<-drop1(mod, scope=
            eval(parse(text=paste("~",Source[term]))))
        df[term]<-mod.1$Df[2]
        SS[term]<-mod.1$Deviance[2] - deviance
        f[term]<-(SS[term]/df[term])/dispersion
        p[term]<-1-pf(f[term], df[term], df.res)
        }
     result<-data.frame(SS, df, f, p)
     row.names(result)<-c(Source, "Residuals")
     names(result)<-c("SS", "Df", "F", "Pr(>F)")
     class(result)<-c("anova","data.frame")
     attr(result,"heading")<-c("Anova Table (Type III tests)\n", paste("Response:", responseName(mod)))
     result
     }
     
        # type II
        
            # Wald test
        
Anova.II.Wald.glm<-function(mod, ...){
    # last modified by J.Fox 4 April 2005
    chisq.term<-function(term){
        which.term<-which(term==names)
        subs.term<-which(assign==which.term)
        relatives<-relatives(term, names, fac)
        subs.relatives<-NULL
        for (relative in relatives) subs.relatives<-c(subs.relatives, which(assign==relative))
        hyp.matrix.1<-I.p[subs.relatives,]
        hyp.matrix.2<-I.p[c(subs.relatives,subs.term),]
        sumry<-summary(mod, corr=FALSE)
        chisq.1<-if (length(subs.relatives)==0) 0 
            else linear.hypothesis(mod, hyp.matrix.1, summary.model=sumry)$Chisq[2]
        chisq.2<-linear.hypothesis(mod, hyp.matrix.2, summary.model=sumry)$Chisq[2]
        chisq.2-chisq.1
        }
    fac<-attr(mod$terms, "factors")
    intercept<-has.intercept(mod)
    p<-length(coefficients(mod))
    I.p<-diag(p)
    names<-term.names(mod)
    if (intercept) names<-names[-1]
    assign<-rep(1:length(names), df.terms(mod))
    assign<-if (intercept) c(0,assign) else assign
    n.terms<-length(names)
    Wald<-rep(0, n.terms)
    df<-rep(0, n.terms)
    p<-rep(0, n.terms)
    for (i in 1:n.terms){
        Wald[i]<-chisq.term(names[i])
        df[i]<-df.terms(mod, names[i])
        p[i]<-1-pchisq(Wald[i],df[i])
        }    
    result<-data.frame(Wald, df, p)
    row.names(result)<-names
    names(result)<-c("Wald Chisq","Df","Pr(>Chisq)")
    class(result)<-c("anova","data.frame")
    attr(result,"heading")<-c("Anova Table (Type II tests)\n", paste("Response:", responseName(mod)))
    result
    }

            # LR test
            
Anova.II.LR.glm <- function(mod, ...){
    # last modified 5 Nov 2002 by J. Fox
    # (some code adapted from drop1.glm)
    which.nms <- function(name) which(asgn == which(names == name))
    fac <- attr(mod$terms, "factors")
    names <- if (has.intercept(mod)) term.names(mod)[-1]
        else term.names(mod)
    n.terms <- length(names)
    X <- model.matrix(mod)
    y <- mod$y
    if (is.null(y)) y <- model.response(model.frame(mod), "numeric")
    wt <- mod$prior.weights
    if (is.null(wt)) wt <- rep(1, length(y))
    asgn <- attr(X, 'assign')
    LR <- rep(0, n.terms)
    df <- df.terms(mod)
    p <- rep(0, n.terms)
    dispersion <- summary(mod, corr = FALSE)$dispersion
    for (term in 1:n.terms){
        rels <- names[relatives(names[term], names, fac)]
        exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
        mod.1 <- glm.fit(X[, -exclude.1, drop = FALSE], y, wt, offset = mod$offset, 
            family = mod$family, control = mod$control)
        dev.1 <- deviance(mod.1)
        mod.2 <- if (length(rels) == 0) mod
            else {
                exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
                glm.fit(X[, -exclude.2, drop = FALSE], y, wt, offset = mod$offset, 
                    family = mod$family, control = mod$control)
                }
        dev.2 <- deviance(mod.2)
        LR[term] <- (dev.1 - dev.2)/dispersion
        p[term] <- 1 - pchisq(LR[term], df[term])
        }
     result <- data.frame(LR, df, p)
     row.names(result) <- names
     names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
     class(result) <- c("anova", "data.frame")
     attr(result,"heading") <- 
        c("Anova Table (Type II tests)\n", paste("Response:", responseName(mod)))
     result
     }


            # F test
            
Anova.II.F.glm <- function(mod, error, error.estimate, ...){
    # last modified 25 Apr 2003 by J. Fox
    # (some code adapted from drop1.glm)
    fam <- family(mod)$family
    if (fam == "binomial" || fam == "poisson") 
        warning("dispersion parameter estimated from the Pearson residuals, not taken as 1")
    which.nms <- function(name) which(asgn == which(names == name))
    if (missing(error)) error <- mod
    df.res <- df.residual(error)
    error.SS <- switch(error.estimate,
        pearson = sum(residuals(error, "pearson")^2),
        dispersion = df.res*summary(error, corr = FALSE)$dispersion,
        deviance = deviance(error))
    fac <- attr(mod$terms, "factors")
    names<-if (has.intercept(mod)) term.names(mod)[-1]
        else term.names(mod)
    n.terms <- length(names)
    X <- model.matrix(mod)
    y <- mod$y
    if (is.null(y)) y <- model.response(model.frame(mod), "numeric")
    wt <- mod$prior.weights
    if (is.null(wt)) wt <- rep(1, length(y))
    asgn <- attr(X, 'assign')
    p <- df <- f <- SS <- rep(0, n.terms+1)
    f[n.terms+1] <- p[n.terms+1] <- NA
    df[n.terms+1] <- df.res
    SS[n.terms+1] <- error.SS
    dispersion <- error.SS/df.res
    df <- c(df.terms(mod), df.res)
    for (term in 1:n.terms){
        rels <- names[relatives(names[term], names, fac)]
        exclude.1 <- as.vector(unlist(sapply(c(names[term], rels), which.nms)))
        mod.1 <- glm.fit(X[, -exclude.1, drop = FALSE], y, wt, offset = mod$offset, 
            family = mod$family, control = mod$control)
        dev.1 <- deviance(mod.1)
        mod.2 <- if (length(rels) == 0) mod
            else {
                exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
                glm.fit(X[, -exclude.2, drop = FALSE], y, wt, offset = mod$offset, 
                    family = mod$family, control = mod$control)
                }
        dev.2 <- deviance(mod.2)
        SS[term] <- dev.1 - dev.2
        f[term] <- SS[term]/(dispersion*df[term])
        p[term] <- 1 - pf(f[term], df[term], df.res)
        }
     result <- data.frame(SS, df, f, p)
     row.names(result) <- c(names, "Residuals")
     names(result) <- c("SS","Df","F","Pr(>F)")
     class(result) <- c("anova","data.frame")
     attr(result,"heading") <- c("Anova Table (Type II tests)\n", 
        paste("Response:", responseName(mod)))
     result
     }

# multinomial logit models (via multinom in the nnet package)
# last modified: 1 Dec 06 by J. Fox

Anova.multinom <-
function (mod, type = c("II", "III", 2, 3), ...)
{
    type <- match.arg(type)
    switch(type,
        II = Anova.II.multinom(mod, ...),
        III = Anova.III.multinom(mod, ...),
        "2" = Anova.II.multinom(mod, ...),
        "3" = Anova.III.multinom(mod, ...))
}

Anova.II.multinom <- function (mod, ...)
{
    which.nms <- function(name) which(asgn == which(names ==
        name))
    fac <- attr(mod$terms, "factors")
    names <- if (has.intercept(mod)) term.names(mod)[-1]
        else term.names(mod)
    n.terms <- length(names)
    X <- model.matrix(mod)
    y <- model.response(model.frame(mod))
    wt <- mod$weights
    asgn <- attr(X, "assign")
    LR <- rep(0, n.terms)
    df <- df.terms(mod)
    p <- rep(0, n.terms)
    for (term in 1:n.terms) {
        rels <- names[relatives(names[term], names, fac)]
        exclude.1 <- as.vector(unlist(sapply(c(names[term], rels),
            which.nms)))
        mod.1 <- multinom(y ~ X[, -c(1, exclude.1)], weights=wt, trace=FALSE)
        dev.1 <- deviance(mod.1)
        mod.2 <- if (length(rels) == 0)
            mod
        else {
            exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
            multinom(y ~ X[, -c(1, exclude.2)], weights=wt, trace=FALSE)
        }
        dev.2 <- deviance(mod.2)
        LR[term] <- dev.1 - dev.2
        p[term] <- 1 - pchisq(LR[term], df[term])
    }
    result <- data.frame(LR, df, p)
    row.names(result) <- names
    names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Anova Table (Type II tests)\n",
        paste("Response:", responseName(mod)))
    result
}

Anova.III.multinom <- function (mod, ...)
{
    names <- if (has.intercept(mod)) term.names(mod)[-1]
        else term.names(mod)
    n.terms <- length(names)
    X <- model.matrix(mod)
    y <- model.response(model.frame(mod))
    wt <- mod$weights
    asgn <- attr(X, "assign")
    LR <- rep(0, n.terms)
    df <- df.terms(mod)
    p <- rep(0, n.terms)
    deviance <- deviance(mod)
    for (term in 1:n.terms) {
        mod.1 <- multinom(y ~ X[, term != asgn][, -1], weights=wt, trace=FALSE)
        LR[term] <- deviance(mod.1) - deviance
        p[term] <- 1 - pchisq(LR[term], df[term])
    }
    result <- data.frame(LR, df, p)
    row.names(result) <- names
    names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Anova Table (Type III tests)\n",
        paste("Response:", responseName(mod)))
    result
}


# proportional-odds logit models (via polr in the MASS package)
# last modified 1 Oct 05 by J. Fox

 Anova.polr <-
function (mod, type = c("II", "III", 2, 3), ...)
{
    type <- match.arg(type)
    switch(type,
        II = Anova.II.polr(mod, ...),
        III = Anova.III.polr(mod, ...),
        "2" = Anova.II.polr(mod, ...),
        "3" = Anova.III.polr(mod, ...))
}

Anova.II.polr <- function (mod, ...)
{
    which.nms <- function(name) which(asgn == which(names ==
        name))
    fac <- attr(mod$terms, "factors")
    names <- term.names(mod)
    n.terms <- length(names)
    X <- model.matrix(mod)
    y <- model.response(model.frame(mod))
    wt <- model.weights(model.frame(mod))
    asgn <- attr(X, "assign")
    LR <- rep(0, n.terms)
    df <- df.terms(mod)
    p <- rep(0, n.terms)
    for (term in 1:n.terms) {
        rels <- names[relatives(names[term], names, fac)]
        exclude.1 <- as.vector(unlist(sapply(c(names[term], rels),
            which.nms)))
        mod.1 <- polr(y ~ X[, -c(1, exclude.1)], weights=wt)
        dev.1 <- deviance(mod.1)
        mod.2 <- if (length(rels) == 0)
            mod
        else {
            exclude.2 <- as.vector(unlist(sapply(rels, which.nms)))
            polr(y ~ X[, -c(1, exclude.2)], weights=wt)
        }
        dev.2 <- deviance(mod.2)
        LR[term] <- dev.1 - dev.2
        p[term] <- 1 - pchisq(LR[term], df[term])
    }
    result <- data.frame(LR, df, p)
    row.names(result) <- names
    names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Anova Table (Type II tests)\n",
        paste("Response:", responseName(mod)))
    result
}

Anova.III.polr <- function (mod, ...)
{
    names <- term.names(mod)
    n.terms <- length(names)
    X <- model.matrix(mod)
    y <- model.response(model.frame(mod))
    wt <- model.weights(model.frame(mod))
    asgn <- attr(X, "assign")
    LR <- rep(0, n.terms)
    df <- df.terms(mod)
    p <- rep(0, n.terms)
    deviance <- deviance(mod)
    for (term in 1:n.terms) {
        mod.1 <- polr(y ~ X[, term != asgn][, -1], weights=wt)
        LR[term] <- deviance(mod.1) - deviance
        p[term] <- 1 - pchisq(LR[term], df[term])
    }
    result <- data.frame(LR, df, p)
    row.names(result) <- names
    names(result) <- c("LR Chisq", "Df", "Pr(>Chisq)")
    class(result) <- c("anova", "data.frame")
    attr(result, "heading") <- c("Anova Table (Type III tests)\n",
        paste("Response:", responseName(mod)))
    result
}

# multivariate linear models
# last modified 1 Dec 06 by J. Fox

has.intercept.mlm <- function (model, ...) 
    any(row.names(coefficients(model)) == "(Intercept)")

Anova.mlm <- function(mod, type=c("II","III", 2, 3), SSPE, error.df, idata, 
    idesign, icontrasts=c("contr.sum", "contr.poly"),
    test.statistic=c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),...){
    type <- match.arg(type)
    test.statistic <- match.arg(test.statistic)
    if (missing(SSPE)) SSPE <- crossprod(residuals(mod))
    if (missing(idata)) {
        idata <- NULL
        idesign <- NULL
        }
    error.df <- if (missing(error.df)) df.residual(mod)
        else error.df
    switch(type,
        II=Anova.II.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, test.statistic, ...),
        III=Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, test.statistic, ...),
        "2"=Anova.II.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, test.statistic, ...),
        "3"=Anova.III.mlm(mod, SSPE, error.df, idata, idesign, icontrasts, test.statistic, ...))
    }
    
Anova.III.mlm <- function(mod, SSPE, error.df, idata, idesign, icontrasts, test, ...){
    intercept <- has.intercept(mod)
    V <- solve(crossprod(model.matrix(mod)))
    p <- nrow(coefficients(mod))
    I.p <- diag(p)
    terms <- term.names(mod)
    n.terms <- length(terms)
    assign <- mod$assign
    if (is.null(idata)){
        if ((n.terms == 0) && intercept) {
             Test <- linear.hypothesis(mod, 1, SSPE=SSPE, ...)
             result <- list(SSP=Test$SSPH, SSPE=SSPE, df=1, error.df=error.df,
                terms="(Intercept)", repeated=FALSE, type="III", test=test)
             class(result) <- "Anova.mlm"
             return(result)
             }
        SSP <- as.list(rep(0, n.terms))
        df <- rep(0, n.terms)
        names(df) <- names(SSP) <- terms
        for (term in 1:n.terms){
            subs <- which(assign == term - intercept)
            hyp.matrix <- I.p[subs,]
            Test <- linear.hypothesis(mod, hyp.matrix, SSPE=SSPE, ...)
            SSP[[term]] <- Test$SSPH
            df[term]<- length(subs)
            }
         result <- list(SSP=SSP, SSPE=SSPE, df=df, error.df=error.df, terms=terms,
            repeated=FALSE, type="III", test=test)
         }
    else {
        if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
        X.design <- model.matrix(idesign, data=idata)
        intercept <- has.intercept(X.design)
        iterms <- term.names(idesign)
        if (intercept) iterms <- c("(Intercept)", iterms)
        df <- rep(0, n.terms*length(iterms))
        hnames <- rep("", length(df))
        P <- SSPEH <- SSP <- as.list(df)
        i <- 0
        for (iterm in iterms){
            for (term in 1:n.terms){
                subs <- which(assign == term - intercept)
                hyp.matrix <- I.p[subs,]
                i <- i + 1
                Test <- linear.hypothesis(mod, hyp.matrix, SSPE=SSPE, 
                    idata=idata, idesign=idesign, icontrasts=icontrasts, iterms=iterm, ...)
                SSP[[i]] <- Test$SSPH
                SSPEH[[i]] <- Test$SSPE
                P[[i]] <- Test$P
                df[i] <- length(subs)
                hnames[i] <- if (iterm == "(Intercept)") terms[term]
                    else if (terms[term] == "(Intercept)") iterm
                    else paste(terms[term], ":", iterm, sep="")
                }
            }
        names(df) <- names(SSP) <- names(SSPEH) <- hnames
        result <- list(SSP=SSP, SSPE=SSPEH, P=P, df=df, error.df=error.df,
            terms=hnames, repeated=TRUE, type="III", test=test)       
        }
     class(result) <- "Anova.mlm"
     result
     }
     
Anova.II.mlm <- function(mod, SSPE, error.df, idata, idesign, icontrasts, test, ...){
    V <- solve(crossprod(model.matrix(mod)))
    SSP.term <- function(term, iterm){
        which.term <- which(term == terms)
        subs.term <- which(assign == which.term)
        relatives <- relatives(term, terms, fac)
        subs.relatives <- NULL
        for (relative in relatives) subs.relatives <- c(subs.relatives, which(assign==relative))
        hyp.matrix.1 <- I.p[subs.relatives,]
        hyp.matrix.2 <- I.p[c(subs.relatives, subs.term),]
        if (missing(iterm)){
            SSP1 <- if (length(subs.relatives) == 0) 0 
                else linear.hypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, ...)$SSPH
            SSP2 <- linear.hypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, ...)$SSPH
            return(SSP2 - SSP1)
            }
        else {
            SSP1 <- if (length(subs.relatives) == 0) 0 
                else linear.hypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
                    idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, ...)$SSPH
            lh2 <- linear.hypothesis(mod, hyp.matrix.2, SSPE=SSPE, V=V, 
                idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, ...)
            return(list(SSP = lh2$SSPH - SSP1, SSPE=lh2$SSPE, P=lh2$P))
            }
        }
    fac <- attr(mod$terms, "factors")
    intercept <- has.intercept(mod)
    p <- nrow(coefficients(mod))
    I.p <- diag(p)
    assign <- mod$assign
    terms <- term.names(mod)
    if (intercept) terms <- terms[-1]
    n.terms <- length(terms)
    if (is.null(idata)){
        if ((n.terms == 0) && intercept) {
             Test <- linear.hypothesis(mod, 1, SSPE=SSPE, ...)
             result <- list(SSP=list(Test$SSPH), SSPE=SSPE, df=1, error.df=error.df,
                terms="(Intercept)", repeated=FALSE, type="II", test=test)
             class(result) <- "Anova.mlm"
             return(result)
             }
        SSP <- as.list(rep(0, n.terms))
        df <- rep(0, n.terms)
        names(df) <- names(SSP) <- terms
        for (i in 1:n.terms){
            SSP[[i]] <- SSP.term(terms[i])
            df[i]<- df.terms(mod, terms[i])
            }    
         result <- list(SSP=SSP, SSPE=SSPE, df=df, error.df=error.df, terms=terms,
            repeated=FALSE, type="II", test=test)
         }
    else {
        if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
        X.design <- model.matrix(idesign, data=idata)
        iintercept <- has.intercept(X.design)
        iterms <- term.names(idesign)
        if (iintercept) iterms <- c("(Intercept)", iterms)
        df <- rep(0, (n.terms + intercept)*length(iterms))
        hnames <- rep("", length(df))
        P <- SSPEH <- SSP <- as.list(df)
        i <- 0
        for (iterm in iterms){
            if (intercept){
                i <- i + 1
                hyp.matrix.1 <- I.p[-1,]
                SSP1 <- linear.hypothesis(mod, hyp.matrix.1, SSPE=SSPE, V=V, 
                    idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, ...)$SSPH
                lh2 <- linear.hypothesis(mod, I.p, SSPE=SSPE, V=V, 
                    idata=idata, idesign=idesign, iterms=iterm, icontrasts=icontrasts, ...)
                SSP[[i]] <- lh2$SSPH - SSP1
                SSPEH[[i]] <- lh2$SSPE
                P[[i]] <- lh2$P
                df[i] <- 1
                hnames[i] <- iterm
                }
            for (term in 1:n.terms){
                subs <- which(assign == term)
                i <- i + 1
                Test <- SSP.term(terms[term], iterm)
                SSP[[i]] <- Test$SSP
                SSPEH[[i]] <- Test$SSPE
                P[[i]] <- Test$P
                df[i]<- length(subs)
                hnames[i] <- if (iterm == "(Intercept)") terms[term]
                    else paste(terms[term], ":", iterm, sep="")
                }
            }
        if (intercept){
            SSP <- SSP[-1]
            SSPEH <- SSPEH[-1]
            P <- P[-1]
            df <- df[-1]
            hnames <- hnames[-1]
            }
        names(df) <- names(P) <- names(SSP) <- names(SSPEH) <- hnames
        result <- list(SSP=SSP, SSPE=SSPEH, P=P, df=df, error.df=error.df,
            terms=hnames, repeated=TRUE, type="II", test=test)       
        }
     class(result) <- "Anova.mlm"
     result
    }

     
print.Anova.mlm <- function(x, ...){
    test <- x$test
    repeated <- x$repeated
    ntests <- length(x$terms)
    tests <- matrix(NA, ntests, 4)
    if (!repeated) SSPE.qr <- qr(x$SSPE) 
    for (term in 1:ntests){
    # some of the code here adapted from stats:::summary.manova
        eigs <- Re(eigen(qr.coef(if (repeated) qr(x$SSPE[[term]]) else SSPE.qr,
            x$SSP[[term]]), symmetric = FALSE)$values)
        tests[term, 1:4] <- switch(test,
            Pillai = stats:::Pillai(eigs, x$df[term], x$error.df),
            Wilks = stats:::Wilks(eigs, x$df[term], x$error.df),
            "Hotelling-Lawley" = stats:::HL(eigs, x$df[term], x$error.df),
            Roy = stats:::Roy(eigs, x$df[term], x$error.df))
        }
    ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
    ok <- !is.na(ok) & ok
    tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3], tests[ok, 4], 
        lower.tail = FALSE))
    rownames(tests) <- x$terms
    colnames(tests) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)")
    tests <- structure(as.data.frame(tests), 
        heading = paste("\nType ", x$type, if (repeated) " Repeated Measures",
            " MANOVA Tests: ", test, " test statistic", sep=""), 
        class = c("anova", "data.frame"))
    print(tests)      
    invisible(x)
    }
    
summary.Anova.mlm <- function(object, test.statistic, multivariate=TRUE, univariate=TRUE, 
    digits=unlist(options("digits")), ...){
    GG <- function(SSPE, P){ # Greenhouse-Geisser correction
        p <- nrow(SSPE)
        if (p < 2) return(NA) 
        lambda <- eigen(SSPE %*% solve(t(P) %*% P))$values
        lambda <- lambda[lambda > 0]
        ((sum(lambda)/p)^2)/(sum(lambda^2)/p)
        }
    HF <- function(gg, error.df, p){ # Huynh-Feldt correction
        ((error.df + 1)*p*gg - 2)/(p*(error.df - p*gg))
        }        
    if (missing(test.statistic)) test.statistic <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
    test.statistic <- match.arg(test.statistic, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
            several.ok=TRUE)
    nterms <- length(object$terms)
    if (multivariate || !object$repeated){       
        cat(paste("\nType ", object$type, if (object$repeated) " Repeated Measures",
                    " MANOVA Tests:\n", sep=""))
        if (!object$repeated){ 
            cat("\nSum of squares and products for error:\n")
            print(object$SSPE, digits=digits)
            }
        for (term in 1:nterms){
            cat(paste("\n------------------------------------------\n",
                "\nTerm:", object$terms[term], "\n"))
            hyp <- list(SSPH=object$SSP[[term]], 
                SSPE=if (object$repeated) object$SSPE[[term]] else object$SSPE,
                P=if (object$repeated) object$P[[term]] else NULL, 
                test=test.statistic, df=object$df[term], 
                df.residual=object$error.df, title=object$terms[term])
            class(hyp) <- "linear.hypothesis.mlm"
            print(hyp, digits=digits, SSPE=object$repeated, ...)
            }
        }
    if (object$repeated && univariate){
        error.df <- object$error.df
        table <- matrix(0, nterms, 6)
        table2 <- matrix(0, nterms, 4)
        rownames(table2) <- rownames(table) <- object$terms
        colnames(table) <- c("SS", "num Df", "Error SS", "den Df", "F", "Pr(>F)")
        colnames(table2) <- c("GG eps", "Pr(>F[GG])",  "HF eps", "Pr(>F[HF])")
        for (term in 1:nterms){
            SSP <- object$SSP[[term]]
            SSPE <- object$SSPE[[term]]
            P <- object$P[[term]]
            p <- ncol(P)
            PtPinv <- solve(t(P) %*% P)
            gg <- GG(SSPE, P)
            table[term, "SS"] <- sum(diag(SSP %*% PtPinv))
            table[term, "Error SS"] <- sum(diag(SSPE %*% PtPinv))
            table[term, "num Df"] <- object$df[term] * p
            table[term, "den Df"] <- error.df * p
            table[term, "F"] <-  (table[term, "SS"]/table[term, "num Df"])/
                (table[term, "Error SS"]/table[term, "den Df"])
            table[term, "Pr(>F)"] <- pf(table[term, "F"], table[term, "num Df"],
                table[term, "den Df"], lower=FALSE)
            table2[term, "GG eps"] <- gg
            table2[term, "HF eps"] <- HF(gg, error.df, p) 
            }
        cat("\nUnivariate Type", object$type, 
            "Repeated-Measures ANOVA Assuming Compound Symmetry\n\n")
        print.anova(table)
        cat("\n\nGreenhouse-Geisser and Huynh-Feldt Corrections\n",
            "for Departure from Compound Symmetry\n\n")
        table2[,"Pr(>F[GG])"] <- pf(table[,"F"], table2[,"GG eps"]*table[,"num Df"],
                table2[,"GG eps"]*table[,"den Df"], lower=FALSE)
        table2[,"Pr(>F[HF])"] <- pf(table[,"F"], table2[,"HF eps"]*table[,"num Df"],
                table2[,"HF eps"]*table[,"den Df"], lower=FALSE)
        table2 <- na.omit(table2)
        print.anova(table2[,1:2])
        cat("\n")
        print.anova(table2[,3:4])
        }
    invisible(object)
    }
    
Anova.manova <- function(mod, ...){
    class(mod) <- c("mlm", "lm")
    Anova(mod, ...)
    }

Manova <- function(mod, ...){
    UseMethod("Manova")
    }
    
Manova.mlm <- function(mod, ...){
    Anova(mod, ...)
    }