# Linear hypothesis tests for lm and glm (J. Fox)
# last modified 9 Nov 02
# 2005-03-17: modified by Achim Zeleis to include `test' and `vcov' arguments,
#             and prettified output


lht<-function(...) linear.hypothesis(...)

linear.hypothesis<-function (model, ...) {
    UseMethod("linear.hypothesis")
    }
    
linear.hypothesis.lm<-function(model, hypothesis.matrix, rhs=0, 
        summary.model=summary(model, corr = FALSE),
	test=c("F", "Chisq"), vcov=NULL,
        white.adjust=FALSE, error.SS, error.df, ...) {
    if (is.aliased(model)) stop ("One or more terms aliased in model.")

    ## obtain RSS and df
    df <- if (missing(error.df)) model$df.residual else error.df
    error.SS <- if(missing(error.SS)) summary.model$sigma^2 * df
    s2 <- error.SS/df

    ## obtain covariance matrix:
    ## can be NULL (default: then white.adjust is used as before),
    ## function or matrix
    if(identical(white.adjust, TRUE)) white.adjust <- "hc3"
    V<-if(is.null(vcov)) {
         if (!is.character(white.adjust)) s2 * summary.model$cov.unscaled
         else hccm(model, type=white.adjust)
       } else {
         if(is.function(vcov)) vcov(model)
	 else vcov
       }

    b <- coefficients(model)
    L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
         else hypothesis.matrix
    q <- nrow(L)

    SSH <- as.vector(t(L %*% b - rhs) %*% inv(L %*% V %*% t(L)) %*% (L %*% b - rhs))

    ## support both asymptotic Chisq test and (approximate) finite sample F test
    ## if df not finite, choose to use Chisq test
    test <- match.arg(test)
    if(!(is.finite(df) && df > 0)) test <- "Chisq"

    ## put together output
    title <- "Linear hypothesis test\n"
    topnote <- paste("Model 1: ", deparse(formula(model)), "\n",
                     "Model 2: restricted model", sep = "")
    rval <- matrix(rep(NA, 8), ncol = 4)
    colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
    rownames(rval) <- 1:2
    rval[,1] <- c(df, df+q)
    
    if(test == "F") {
      f <- SSH/q
      p <- pf(f, q, df, lower.tail = FALSE)
      rval[2,2:4] <- c(-q, f, p)
    } else {
      p <- pchisq(SSH, q, lower.tail = FALSE)
      rval[2,2:4] <- c(-q, SSH, p)
    }
    
    ## if default vcov estimate was used, add RSS to output
    if(is.null(vcov)) {
      rval2 <- matrix(rep(NA, 4), ncol = 2)
      colnames(rval2) <- c("RSS", "Sum of Sq")
      rval2[,1] <- c(error.SS, error.SS + SSH * s2)
      rval2[2,2] <- -diff(rval2[,1])
      rval <- cbind(rval, rval2)[,c(1, 5, 2, 6, 3, 4)]
    }    
    
    rval <- structure(as.data.frame(rval), heading = c(title, topnote),
	    class = c("anova", "data.frame"))
    return(rval)
    }

linear.hypothesis.glm<-function(model, hypothesis.matrix, rhs=0, 
        summary.model=summary(model, corr = FALSE),
	test=c("Chisq", "F"), vcov=NULL, error.df, ...) {
    if (is.aliased(model)) stop ("One or more terms aliased in model.")

    ## obtain df
    df <- if (missing(error.df)) model$df.residual else error.df

    ## compute covariance matrix (see linear.hypothesis.lm)
    V<-if(is.null(vcov)) {
         summary.model$dispersion*summary.model$cov.unscaled
       } else {
         if(is.function(vcov)) vcov(model)
	 else vcov
       }
    
    b <- coefficients(model)
    L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
         else hypothesis.matrix
    q <- nrow(L)
    SSH <- as.vector(t(L %*% b - rhs) %*% inv(L %*% V %*% t(L)) %*% (L %*% b - rhs))

    ## support both asymptotic Chisq test and (approximate) finite sample F test
    ## if not finite, choose to use Chisq test
    test <- match.arg(test)
    if(!(is.finite(df) && df > 0)) test <- "Chisq"

    ## set up return value
    title <- "Linear hypothesis test\n"
    topnote <- paste("Model 1: ", deparse(formula(model)), "\n",
                     "Model 2: restricted model", sep = "")
    rval <- matrix(rep(NA, 8), ncol = 4)
    colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
    rownames(rval) <- 1:2
    rval[,1] <- c(df, df+q)
    
    if(test == "F") {
      f <- SSH/q
      p <- pf(f, q, df, lower.tail = FALSE)
      rval[2,2:4] <- c(-q, f, p)
    } else {
      p <- pchisq(SSH, q, lower.tail = FALSE)
      rval[2,2:4] <- c(-q, SSH, p)
    }

    rval <- structure(as.data.frame(rval), heading = c(title, topnote),
	    class = c("anova", "data.frame"))
    return(rval)
    }
