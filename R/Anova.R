# Type II and III tests for linear and generalized linear models (J. Fox)

# last modified 31 Jan 04

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
 
Anova.lm<-function(mod, error, type=c("II","III"), ...){
    type<-match.arg(type)
    switch(type,
        II=Anova.II.lm(mod, error, ...),
        III=Anova.III.lm(mod, error, ...))
    }

Anova.aov <- function(mod, ...){
    # last modified 8 Mar 2002 by J. Fox
    class(mod) <- "lm"
    Anova.lm(mod, ...)
    }

        # type II
        
Anova.II.lm<-function(mod, error, ...){
    # last modified by J.Fox 11 Dec 2000
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
            else linear.hypothesis(mod, hyp.matrix.1, summary.model=sumry, ...)$SSH
        SS2<-linear.hypothesis(mod, hyp.matrix.2, summary.model=sumry, ...)$SSH
        SS2-SS1
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
    # last modified by J.Fox 11 Dec 2000
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
        SS[term]<-test$SSH
        df[term]<-test$Df[1]
        f[term]<-test$f
        p[term]<-test$p
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
    
Anova.glm<-function(mod, type=c("II","III"), test.statistic=c("LR", "Wald", "F"), 
    error, error.estimate=c("pearson", "dispersion", "deviance"), ...){
    #last modified by J. Fox 15 Feb 2001
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
        Wald[term]<-test$ChiSquare
        df[term]<-test$Df
        p[term]<-test$p
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
    # last modified by J.Fox 11 Dec 2000
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
            else linear.hypothesis(mod, hyp.matrix.1, summary.model=sumry)$ChiSquare
        chisq.2<-linear.hypothesis(mod, hyp.matrix.2, summary.model=sumry)$ChiSquare
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
