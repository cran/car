# score test of nonconstant variance (J. Fox)

ncv.test<-function(model, ...){
    # last modified 15 Dec 2000 by J. Fox
    UseMethod("ncv.test")
    }

ncv.test.lm<-function (model, var.formula, data=NULL, subset, na.action) {
    # last modified 29 July 2001 by J. Fox
    if (!is.null(weights(model))) stop("requires unweighted linear model")
    sumry<-summary(model)
    residuals<-residuals(model)
    S.sq<-df.residual(model)*(sumry$sigma)^2/sum(!is.na(residuals))
    U<-(residuals^2)/S.sq
    if (missing(var.formula)) {
        mod<-lm(U~fitted.values(model))
        varnames<-"fitted.values"
        var.formula<-~fitted.values
        df<-1
        }
    else {
        if (missing(na.action)) 
            na.action <- options()$na.action
        m <- match.call(expand.dots = FALSE)
        if (is.matrix(eval(m$data, sys.frame(sys.parent())))) 
            m$data <- as.data.frame(data)
        m$formula<-var.formula
        m$var.formula <- m$model <- m$... <- NULL
        m[[1]] <- as.name("model.frame")
        mf <- eval(m, sys.frame(sys.parent()))
        response <- attr(attr(mf, "terms"), "response")
        if (response) stop(paste("Variance formula contains a response."))
        mf$U<-U
        .X<-model.matrix(as.formula(paste("U~",as.character(var.formula)[2],"-1")), data=mf)
        mod<-lm(U~.X)
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
