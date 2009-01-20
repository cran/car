# score test of nonconstant variance (J. Fox)

# last modified 20 January 2009 by J. Fox

ncv.test<-function(model, ...){
    # last modified 15 Dec 2000 by J. Fox
    UseMethod("ncv.test")
    }

ncv.test.lm<-function (model, var.formula, data=NULL, subset, na.action, ...) {
    # last modified 16 Jan 2005 by J. Fox
    if ((!is.null(class(model$na.action))) && class(model$na.action) == 'exclude') 
        model <- update(model, na.action=na.omit)
    sumry<-summary(model)
    residuals <- residuals(model, type="pearson") # suggested by S. Weisberg
    S.sq<-df.residual(model)*(sumry$sigma)^2/sum(!is.na(residuals))
    U<-(residuals^2)/S.sq
    if (missing(var.formula)) {
        mod<-lm(U~fitted.values(model))
        varnames<-"fitted.values"
        var.formula<-~fitted.values
        df<-1
        }
    else {
        if (missing(na.action)){
            na.action <- if (is.null(model$na.action)) options()$na.action
                else parse(text=paste('na.',class(model$na.action), sep=''))
            }
        m <- match.call(expand.dots = FALSE)
        if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
            m$data <- as.data.frame(data)
        m$formula<-var.formula
        m$var.formula <- m$model <- m$... <- NULL
        m[[1]] <- as.name("model.frame")
        mf <- eval(m, sys.frame(sys.parent()))
        response <- attr(attr(mf, "terms"), "response")
        if (response) stop(paste("Variance formula contains a response."))
        .X<-model.matrix(as.formula(paste("~",as.character(var.formula)[2],"-1")), data=mf)
        common.obs <- intersect(names(U), rownames(.X))
        mod<-lm(U[common.obs]~.X[common.obs,])
        df<-sum(!is.na(coefficients(mod)))-1
        }
    SS<-anova(mod)$"Sum Sq"
    RegSS<-sum(SS)-SS[length(SS)]
    Chisq<-RegSS/2
    result<-list(formula=var.formula, formula.name="Variance", ChiSquare=Chisq, Df=df, 
        p=1-pchisq(Chisq, df), test="Non-constant Variance Score Test")
    class(result)<-"chisq.test"
    result
    }
    
ncv.test.glm<-function(model, ...){
    # last modified 15 Dec 2000 by J. Fox
    stop("requires lm object")
    }
    
 print.chisq.test<-function(x, ...){
    title<-if (!is.null(x$test)) x$test else "Chisquare Test"
    cat(title,"\n")
    if (!is.null(x$formula)) cat(x$formula.name, 
        "formula:", as.character(x$formula), "\n")
    cat("Chisquare =", x$ChiSquare,"   Df =", x$Df,
        "    p =", x$p, "\n")
    invisible(x)
    }
