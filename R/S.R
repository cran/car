# Alternatives to `stats` functions for various reasons 12/26/2017
# Summarize.lm:  Adds new argument vcov.=vcov to specify a covariance matrix.  The default reproduces the
#              current output.  The linearHypothesis function is used to compute the overall F-test.
# print.Summarize.lm: new arguments:
#              header=TRUE prints or suppresses the header
#              resid.summary=TRUE prints or suppresses the residual summary
#              adj.r.squared=TRUE prints or suppresses printing of the adjusted r.squared
#              brief=FALSE if TRUE sets header=resid.summary=adj.r.squared=FALSE
#     In addition output is modified to include the vcov. argument if it is not set to vcov
# Confint.lm:  new argument vcov.=vcov where vcov. is either a matrix of the right size or
#     a fuction so that vcov.(object) returns an estmated covariance matrix.
# 2016-12-27 For now, the override functions start with a Capital letter
# 2017-02-10: Renamed using uc letters; introduced default methods. J. Fox
# 2017-02-21: removed Vcov as it is not needed.  Removed vcov=Boot and added an example with
#               b1 <- Boot(object)
#               Summarize(object vcov. = cov(b1$t))
#               Confint(b1)  # to get the same bootstrap and use bca method

# Summarize adds vcov. argument
# 2017-02-23:  S. Weisberg added Summarize.glm and print.Summarize.glm
# 2017-05-15:  S. Weisberg added singular.ook=TRUE to call to linearHypothesis
# 2017-06-15:  S. Weisberg moved arguments from print.Summarize to Summarize
# 2017-06-22:  S. Weisberg added a 'Summarise' method that is the same as 'Summarize'
# 2017-09-20:  J. Fox added estimate and exponentiate arguments to Confint()
# 2017-10-03:  J. Fox fixed bug in Confint.default(), which didn't return its result
#                     added  Confint.polr(), Confint.multinom(), Summarize.multinom(),
#                            print.Summarize.multinom()
# 2017-10-04:  J. Fox added S() generic and methods & tweaked some Summarize() and print() methods
# 2017-10-10:  S. Weisberg fixed bug in dispersion arg in Summarize.glm
# 2017-10-11:  J. Fox modified Confint.glm() to suppress message about profiling likelihood
# 2017-10-12:  J. Fox made changes to Confint.glm() et al. to handle vcov. and dispersion args consistently
# 2017-10-25:  J. Fox added terms and intercept args to S() and methods to print coefficients selectively
# 2017-11-02:  J. Fox added Summarize() methods for lme, lmer, and glmer objects
# 2017-11-07,09:  J. Fox added complete=FALSE to vcov.() calls
# 2017-11-07:  J. Fox added unexported formatCall() for improved formatting of calls
# 2017-11-24:  J. Fox made small improvements to output messages, etc.
# 2017-11-29:  J. Fox made fixes for vcov() and vcov.() calls.
# 2017-12-27:  J. Fox tweaked the Summarize() output for mixed models.
# 2017-12-29:  J. Fox added fit statistics to Summarize() output for various models.
# 2018-01-15:  S. Weisberg all Summmarize/Summarise methods renamed S
# 2018-02-02:  J. Fox fixed S.lm() and S.glm() output when vcov. arg not given.
# 2018-02-07,08,12:  J. Fox removed leading blank lines in formatCall() and elsewhere.
# 2018-10-23: J. Fox made coefs2use() work with models without an intercept even if intercept arg is TRUE.
# 2019-05-02: J. Fox fixed bug in Confint.polr() that exponentiated coefficients twice (reported by Thamron Keowmani).
# 2019-05-02,13: J. Fox made several S() methods tolerant of model with 1 coefficient or
#             in the case of multinom models, 2 response levels(reported by Thamron Keowmani).

formatCall <- function(call){
  call <- if (is.character(call)){
    if (length(call) > 1) paste(call, collapse=" ") else call
  }
  else paste(deparse(call), sep = "", collapse = "")
  call <-  gsub("\\s+", " ", call)
  call <- paste("Call:", call)
  call <- strwrap(call, width=getOption("width"))
  paren <- regexpr("\\(", call[1])
  if (paren > 0 && length(call) > 1){
    call[-1] <- paste0(paste(rep(" ", paren), collapse=""), call[-1])
  }
  paste0(paste(call, collapse="\n"), "\n")
}

fitstats <- function(model){
    logLik <- logLik(model)
    result <- c(logLik=as.vector(logLik), df=attr(logLik, "df"), AIC=AIC(model), BIC=BIC(model))
    class(result) <- "fitstats"
    result
}

print.fitstats <- function(x, digits=2, ...){
    x <- round(x, digits=digits)
    result <- format(x)
    result["df"] <- format(x["df"])
    cat("\n")
    print(result, quote=FALSE)
    invisible(x)
}

S <- function(object, brief, ...){
    UseMethod("S")
}

#Summarise <- function(object, brief, ...){
#    UseMethod("S")
#}

S.default <- function(object, brief, ...) summary(object, ...)

#S.glm <- function(object, ...) {
#  if(object$family$family == "gaussian" & object$family$link == "identity")
#    S.lm(object, ...) else summary(object, ...)
#}

S.lm <- function (object, brief=FALSE, correlation = FALSE, symbolic.cor = FALSE,
                          vcov. = vcov(object, complete=FALSE), header=TRUE, resid.summary=FALSE,
                          adj.r2=FALSE, ...) {
    z <- object
    p <- z$rank
    rdf <- z$df.residual
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        }
        else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        resvar <- rss/rdf
        ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
        class(ans) <- "S.lm"
        ans$aliased <- is.na(coef(object))
        ans$residuals <- r
        ans$df <- c(0L, n, length(ans$aliased))
        ans$coefficients <- matrix(NA, 0L, 4L)
        dimnames(ans$coefficients) <- list(NULL, c("Estimate",
                                                   "Std. Error", "t value", "Pr(>|t|)"))
        ans$sigma <- sqrt(resvar)

        ans$r.squared <- ans$adj.r.squared <- 0
        ans$header <- header
        ans$resid.summary <- resid.summary
        ans$adj.r2 <- adj.r2
        ans$brief <- brief
        ans$fitstats <- round(c(AIC=AIC(object), BIC=BIC(object)), digits=2)
        return(ans)
    }
    if (is.null(z$terms))
        stop("invalid 'lm' object:  no 'terms' component")
    if (!inherits(object, "lm"))
        warning("calling summary.lm(<fake-lm-object>) ...")
    Qr <- object$qr
    n <- NROW(Qr$qr)
    if (is.na(z$df.residual) || n - p != z$df.residual)
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    r <- z$residuals
    f <- z$fitted.values
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept"))
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    if (is.finite(resvar) && resvar < (mean(f)^2 + var(f)) *
        1e-30)
        warning("essentially perfect fit: summary may be unreliable")
    p1 <- 1L:p
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    #    se <- sqrt(diag(R) * resvar)
    V <- if(is.matrix(vcov.)) vcov. else
        if(deparse(substitute(vcov.) == "Boot")) cov((b1 <- Boot(object))$t) else
            vcov.(object)
    se <- sqrt(diag(V))
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval),
                                                    rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]],
                                       c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$aliased <- is.na(coef(object))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
        df.int <- if (attr(z$terms, "intercept"))
            1L
        else 0L
        ans$r.squared <- mss/(mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
        #      ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
        #                          numdf = p - df.int, dendf = rdf)
        # linearHypothesis computes overall F test allowing for alternative covariance matrices
        mat <- diag(p - df.int)
        if(df.int==1) mat <- cbind(0, mat)
        lh <- linearHypothesis(z, mat, vcov.=V, singular.ok=TRUE)
        ans$fstatistic <- c(value = lh$F[2], numdf = lh$Df[2], dendf = lh$Res.Df[2])
    }
    else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1, 1)]
    if (correlation) {
        ans$correlation <- (R * resvar)/outer(se, se)
        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
        ans$symbolic.cor <- symbolic.cor
    }
    if (!is.null(z$na.action))
        ans$na.action <- z$na.action
    ans$vcov. <- if (missing(vcov.)) "" else deparse(substitute(vcov.))
    ans$header <- header
    ans$resid.summary <- resid.summary
    ans$adj.r2 <- adj.r2
    ans$brief <- brief
    ans$fitstats <- round(c(AIC=AIC(object), BIC=BIC(object)), digits=2)
    class(ans) <- "S.lm"
    ans
}

print.S.lm <- function(x, digits = max(3, getOption("digits") - 3),
                               symbolic.cor = x$symbolic.cor,
                               signif.stars = getOption("show.signif.stars"), ...) {
    header <- x$header
    resid.summary <- x$resid.summary
    adj.r2 <- x$adj.r2
    brief <- x$brief
    if (brief) header <- resid.summary <- adj.r2 <- FALSE
    if (header) {
        cat(formatCall(x$call))
        if(x$vcov. != ""){
            cat("Standard errors computed by",  x$vcov., "\n")
        }
    }
    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    if (resid.summary) {
        cat('\n', if (!is.null(x$weights) && diff(range(x$weights)))
            "Weighted ", "Residuals:\n", sep = "")
        if (rdf > 5L) {
            nam <- c("Min", "1Q", "Median", "3Q", "Max")
            rq <- if (length(dim(resid)) == 2L)
                structure(apply(t(resid), 1L, quantile), dimnames = list(nam,
                                                                         dimnames(resid)[[2L]]))
            else {
                zz <- zapsmall(quantile(resid), digits + 1)
                structure(zz, names = nam)
            }
            print(rq, digits = digits, ...)
        }
        else if (rdf > 0L) {
            print(resid, digits = digits, ...)
        }
        else {
            cat("ALL", df[1L], "residuals are 0: no residual degrees of freedom!\n")
        }
    }
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    }
    else {
        if (header || resid.summary) cat("\n")
        if (nsingular <- df[3L] - df[1L])
            cat("Coefficients: (", nsingular, " not defined because of singularities)\n",
                sep = "")
        else cat("Coefficients:\n")
        coefs <- x$coefficients
        if (!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
                                                                    colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                     na.print = "NA", ...)
    }
    cat("\nResidual standard deviation:", format(signif(x$sigma,
                                                        digits)), "on", rdf, "degrees of freedom\n")
    if (nzchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")
    if (!is.null(x$fstatistic)) {
        cat("Multiple R-squared:", formatC(x$r.squared, digits = digits))
        if (adj.r2) {
            cat(",\tAdjusted R-squared:", formatC(x$adj.r.squared,
                                                  digits = digits))
        }
        cat("\nF-statistic:", formatC(x$fstatistic[1L], digits = digits),
            "on", x$fstatistic[2L], "and", x$fstatistic[3L],
            "DF,  p-value:", format.pval(pf(x$fstatistic[1L],
                                            x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE),
                                         digits = digits), "\n")
    }
    print(x$fitstats)
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1L) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.colnames = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2,
                                 digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    cat("\n")
    invisible(x)
}

S.glm <-
    function (object, brief=FALSE, exponentiate, dispersion, correlation = FALSE, symbolic.cor = FALSE,
              vcov. = vcov(object, complete=FALSE), header=TRUE, resid.summary=FALSE, ...)
    {
        vcov.arg <- if (missing(vcov.)) "" else deparse(substitute(vcov.))
        if (missing(exponentiate)) exponentiate <- object$family$link %in% c("log", "logit")
        #    if(!is.null(dispersion)) vcov. <- "vcov" # ignore vcov. arg if dispersion is set
        if (!missing(dispersion) && !missing(vcov.))
            stop("cannot specify both the dispersion and vcov. arguments")
        profile.likelihood <- missing(vcov.) && missing(dispersion)
        est.disp <- FALSE
        df.r <- object$df.residual
        if (missing(dispersion))
            dispersion <- if (object$family$family %in% c("poisson",
                                                          "binomial"))
                1
        else if (df.r > 0) {
            est.disp <- TRUE
            if (any(object$weights == 0))
                warning("observations with zero weight not used for calculating dispersion")
            sum((object$weights * object$residuals^2)[object$weights >
                                                          0])/df.r
        }
        else {
            est.disp <- TRUE
            NaN
        }
        aliased <- is.na(coef(object))
        p <- object$rank
        if (p > 0) {
            p1 <- 1L:p
            Qr <- object$qr
            coef.p <- object$coefficients[Qr$pivot[p1]]
            covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
            dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
            # changed
            covmat <- if(is.matrix(vcov.)) vcov. else
            {if(!est.disp) dispersion * covmat.unscaled else vcov.(object)}
            # end change
            var.cf <- diag(covmat)
            s.err <- sqrt(var.cf)
            tvalue <- coef.p/s.err
            dn <- c("Estimate", "Std. Error")
            if (!est.disp) {
                pvalue <- 2 * pnorm(-abs(tvalue))
                coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
                dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                              "z value", "Pr(>|z|)"))
            }
            else if (df.r > 0) {
                pvalue <- 2 * pt(-abs(tvalue), df.r)
                coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
                dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                              "t value", "Pr(>|t|)"))
            }
            else {
                coef.table <- cbind(coef.p, NaN, NaN, NaN)
                dimnames(coef.table) <- list(names(coef.p), c(dn,
                                                              "t value", "Pr(>|t|)"))
            }
            df.f <- NCOL(Qr$qr)
        }
        else {
            coef.table <- matrix(, 0L, 4L)
            dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
                                                 "t value", "Pr(>|t|)"))
            covmat.unscaled <- covmat <- matrix(, 0L, 0L)
            df.f <- length(aliased)
        }
        keep <- match(c("call", "terms", "family", "deviance", "aic",
                        "contrasts", "df.residual", "null.deviance", "df.null",
                        "iter", "na.action"), names(object), 0L)
        ans <- c(object[keep], list(deviance.resid = residuals(object,
                                                               type = "deviance"), coefficients = coef.table, aliased = aliased,
                                    dispersion = dispersion, df = c(object$rank, df.r, df.f),
                                    cov.unscaled = covmat.unscaled, cov.scaled = covmat))
        if (correlation && p > 0) {
            dd <- sqrt(diag(covmat.unscaled))
            ans$correlation <- covmat.unscaled/outer(dd, dd)
            ans$symbolic.cor <- symbolic.cor
        }
        # add to value
        ans$fitstats <- fitstats(object)
        ans$vcov. <- vcov.arg
        ans$header <- header
        ans$resid.summary <- resid.summary
        ans$brief <- brief
        if (exponentiate) ans$exponentiated <- if (profile.likelihood) Confint(object, exponentiate=TRUE, silent=TRUE)
        else  Confint(object, exponentiate=TRUE, silent=TRUE, vcov.=covmat)
        # end add
        class(ans) <- "S.glm"
        return(ans)
    }

print.S.glm <-
    function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
              signif.stars = getOption("show.signif.stars"), ...)
    {
        header <- x$header
        resid.summary <- x$resid.summary
        brief <- x$brief
        if (brief) {
          header <- resid.summary <-  FALSE
          x$exponentiated <- NULL
        }
        if (header) {
            cat(formatCall(x$call))
            if(x$vcov. != ""){
              cat("Standard errors computed by",  x$vcov., "\n")
            }
        }
        if(resid.summary){
            cat("Deviance Residuals: \n")
            if (x$df.residual > 5) {
                x$deviance.resid <- setNames(quantile(x$deviance.resid,
                                                      na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
            }
            xx <- zapsmall(x$deviance.resid, digits + 1L)
            print.default(xx, digits = digits, na.print = "", print.gap = 2L)
        }
        if (length(x$aliased) == 0L) {
            cat("\nNo Coefficients\n")
        }
        else {
            if (header || resid.summary) cat("\n")
            df <- if ("df" %in% names(x))
                x[["df"]]
            else NULL
            if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
                cat("Coefficients: (", nsingular, " not defined because of singularities)\n",
                    sep = "")
            else cat("Coefficients:\n")
            coefs <- x$coefficients
            if (!is.null(aliased <- x$aliased) && any(aliased)) {
                cn <- names(aliased)
                coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn,
                                                                         colnames(coefs)))
                coefs[!aliased, ] <- x$coefficients
            }
            printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                         na.print = "NA", ...)
        }
        cat("\n(Dispersion parameter for ", x$family$family, " family taken to be ",
            format(x$dispersion), ")\n\n", apply(cbind(paste(format(c("Null",
                                                                      "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance",
                                                                                                                                       "deviance")]), digits = max(5L, digits + 1L)), " on",
                                                       format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"),
                                                 1L, paste, collapse = " "), sep = "")
        if (nzchar(mess <- naprint(x$na.action)))
            cat("  (", mess, ")\n", sep = "")
        print(x$fitstats)
#        cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)), "\n\n",
         cat("\nNumber of Fisher Scoring iterations: ", x$iter,
            "\n", sep = "")
        correl <- x$correlation
        if (!is.null(correl)) {
            p <- NCOL(correl)
            if (p > 1) {
                cat("\nCorrelation of Coefficients:\n")
                if (is.logical(symbolic.cor) && symbolic.cor) {
                    print(symnum(correl, abbr.colnames = NULL))
                }
                else {
                    correl <- format(round(correl, 2L), nsmall = 2L,
                                     digits = digits)
                    correl[!lower.tri(correl)] <- ""
                    print(correl[-1, -p, drop = FALSE], quote = FALSE)
                }
            }
        }
        cat("\n")
        if (!is.null(x$exponentiated)){
            cat("Exponentiated Coefficients and Confidence Bounds\n")
            print(x$exponentiated)
            cat("\n")
        }
        invisible(x)
    }

S.multinom <- function(object, brief=FALSE, exponentiate=FALSE, ...){
    result <- summary(object, ...)
    result$brief <- brief
    result$fitstats <- fitstats(object)
    if (exponentiate) result$exponentiated <- Confint(object, exponentiate=TRUE)
    class(result) <- "S.multinom"
    result
}

print.S.multinom <- function (x, digits = max(3, getOption("digits") - 3),
                                      signif.stars = getOption("show.signif.stars"), ...) {
    if (!x$brief) cat(formatCall(x$call))
    cat("\nCoefficients:\n")
    b <- x$coefficients
    se <- x$standard.errors
    z <- b/se
    p <- 2*pnorm(abs(z), lower.tail=FALSE)
    levels <- x$lev
    if (length(levels) == 2){
        table <- cbind(b, se, z, p)
        colnames(table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
        cat("\n ", levels[2], "\n")
        printCoefmat(table, signif.stars=signif.stars, digits=digits, ...)
    }
    else{
        table <-  abind(t(b), t(se), t(z), t(p), along=1.5)
        dimnames(table)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
        for (level in levels[-1]){
            cat("\n ", level, "\n")
            tab <- table[, , level]
            if (is.vector(tab)){
              cnames <- names(tab)
              tab <- matrix(tab, nrow=1)
              colnames(tab) <- cnames
              rownames(tab) <- x$coefnames
            }
            printCoefmat(tab, signif.stars=signif.stars, digits=digits, ...)
        }
    }
    cat("\nResidual Deviance:", format(x$deviance, digits=digits, ...), "\n")
    print(x$fitstats)
    exponentiated <- x$exponentiated
    if (!is.null(exponentiated)){
        cat("\nExponentiated Coefficients:\n")
      if (length(dim(table)) == 2) print(exponentiated, digits=digits, ...)
      else  for (response in dimnames(table)[[3]]){
            cat("\n ", response, "\n")
            print(exponentiated[, , response], digits=digits, ...)
        }
    }
    invisible(x)
}

S.polr <- function(object, brief=FALSE, exponentiate=FALSE, ...){
    sumry <- summary(object, ...)
    sumry$brief <- brief
    sumry$fitstats <- fitstats(object)
    if (exponentiate){
        sumry$exponentiated <- Confint(object, exponentiate=TRUE, ...)
    }
    class(sumry) <- c("S.polr", class(sumry))
    sumry
}

print.S.polr <- function(x, digits = max(3, getOption("digits") - 3),
                                 signif.stars = getOption("show.signif.stars"), ...) {
    if (!x$brief) cat(formatCall(x$call))
    table <- x$coefficients
    table <- cbind(table, 2*pnorm(abs(table[, 3]), lower.tail=FALSE))
    n.par <- nrow(table)
    n.ints <- length(x$zeta)
    n.coefs <- n.par - n.ints
    coef.table <- table[1:n.coefs, , drop=FALSE]
    int.table <- table[(n.coefs + 1):n.par, , drop=FALSE]
    colnames(coef.table) <- colnames(int.table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    if (!x$brief) cat("\n")
    cat(" Coefficients:\n")
    printCoefmat(coef.table, digits=digits, signif.stars=signif.stars, ...)
    cat("\n Intercepts (Thresholds):\n")
    printCoefmat(int.table, digits=digits, signif.stars=signif.stars, ...)
    cat("\nResidual Deviance:", format(x$deviance), "\n")
    print(x$fitstats)
    if (!is.null(x$exponentiated)){
        cat("\n Exponentiated Coefficients\n")
        print(x$exponentiated)
    }
    invisible(x)
}

S.lmerMod <- function(object, brief=FALSE, KR=FALSE, correlation=FALSE, ...){
  sumry <- summary(object)
  coefs <- sumry$coefficients
  REML <- grepl("REML", sumry$methTitle)
  # the following code for no of groups and obs borrowed from print.merMod()
  dims <- object@devcomp$dims
  ngrps <- vapply(object@flist, nlevels, 0L)
  if (KR){
    if (!REML) stop("KR tests available only for REML estimates")
    b <- coefs[, 1]
    vcov <- as.matrix(pbkrtest::vcovAdj(object))
    coefs[, 2] <- sqrt(diag(vcov))
    p <- length(b)
    coefs <- cbind(coefs, matrix(0, p, 2))
    I.p <- diag(p)
    for (i in 1:p){
      test <- pbkrtest::KRmodcomp(object, I.p[i, , drop=FALSE])
      coefs[i, 3] <- sign(coefs[i, 1])*sqrt(pbkrtest::getKR(test, "Fstat"))
      coefs[i, 4] <- pbkrtest::getKR(test, "ddf")
      coefs[i, 5] <- pbkrtest::getKR(test, "p.value")
    }
    colnames(coefs) <- c("Estimate", "Std. Error", "t value", "df for t", "Pr(>|t|)")
  }
  else {
    coefs <- cbind(coefs, 2*pnorm(abs(coefs[, 3]), lower.tail=FALSE))
    colnames(coefs)[3:4] <- c("z value", "Pr(>|z|)")
    vcov <- as.matrix(sumry$vcov)
  }
  result <- list(logLik=sumry$logLik, fixed.effects=coefs, random.effects=sumry$varcor,
                 REML=REML, KR=KR, call=sumry$call, brief=brief,
                 vcov=vcov, correlation=correlation, nobs=dims[["n"]], ngrps=ngrps,
                 fitstats=fitstats(object))
  class(result) <- "S.lmerMod"
  result
}

print.S.lmerMod <- function(x, digits=max(3, getOption("digits") - 3),
                                 signif.stars = getOption("show.signif.stars"), ...){
  if (!x$brief) {
      cat(paste("Linear mixed model fit by", if (x$REML) "REML" else "ML", "\n"))
      cat(formatCall(x$call))
  }
  if (x$KR) cat("\nEstimates of Fixed Effects with KR Tests\n")
  else cat("\nEstimates of Fixed Effects:\n")
  printCoefmat(x$fixed.effects, digits=digits, signif.stars=signif.stars)
  if (x$correlation) {
    # the following code adapted from print.summary.merMod()
    cor <- cov2cor(x$vcov)
    p <- ncol(cor)
    if (p > 1) {
      rn <- rownames(x$fixed.effects)
      rns <- abbreviate(rn, minlength = 11)
      cat("\nCorrelations of Fixed Effects:\n")
      cor <- matrix(format(round(cor, 3), nsmall = 3),
                    ncol = p, dimnames = list(rns, abbreviate(rn, minlength = 6)))
      cor[!lower.tri(cor)] <- ""
      print(cor[-1, -p, drop = FALSE], quote = FALSE)
    }
  }
  cat("\nEstimates of Random Effects (Covariance Components):\n")
  print(x$random.effects, digits=digits)
  # cat(paste0("\nLog-likelihood (", if (x$REML) "REML) = " else "ML) = ",
  #            format(x$logLik, digits=digits), "\n"))
  # the following code adapted from lme4:::.prt.grps()
  cat(sprintf("\nNumber of obs: %d, groups: ", x$nobs),
      paste(paste(names(x$ngrps), x$ngrps, sep = ", "), collapse = "; "), fill = TRUE)
  print(x$fitstats)
  invisible(x)
}

S.lme <- function(object, brief=FALSE, correlation=FALSE, ...){
  sumry <- summary(object)
  coefs <- sumry$tTable
  colnames(coefs) <- c("Estimate", "Std.Error", "df", "t value", "Pr(>|t|)")
  REML <- sumry$method == "REML"
  result <- list(logLik=sumry$logLik, fixed.effects=coefs,
                 random.effects=summary(sumry$modelStruct),
                 REML=REML, call=sumry$call, brief=brief,
                 vcov=sumry$varFix, sigma=sumry$sigma, dims=sumry$dims,
                 correlation=correlation, fitstats=fitstats(object))
  class(result) <- "S.lme"
  result
}

print.S.lme <- function(x, digits=max(3, getOption("digits") - 3),
                                signif.stars = getOption("show.signif.stars"), ...){
  if (!x$brief) {
      cat(paste("Linear mixed model fit by", if (x$REML) "REML" else "ML"))
      if (!is.null(x$call$data)) cat(",  Data:", as.character(x$call$data))
      cat("\n")
  }

  cat("\nFixed Effects:\n")
  if (!x$brief) cat(" Formula:", deparse(x$call$fixed), "\n\n")
  printCoefmat(x$fixed.effects, digits=digits, signif.stars=signif.stars, cs.ind=1:2)
  if (x$correlation) {
    # the following code adapted from print.summary.merMod()
    cor <- cov2cor(x$vcov)
    p <- ncol(cor)
    if (p > 1) {
      rn <- rownames(x$fixed.effects)
      rns <- abbreviate(rn, minlength = 11)
      cat("\nCorrelations of Fixed Effects:\n")
      cor <- matrix(format(round(cor, 3), nsmall = 3),
                    ncol = p, dimnames = list(rns, abbreviate(rn, minlength = 6)))
      cor[!lower.tri(cor)] <- ""
      print(cor[-1, -p, drop = FALSE], quote = FALSE)
    }
  }
  cat("\n")
  print(x$random.effects, sigma=x$sigma, digits=digits)
  # cat(paste0("\nLog-likelihood (", if (x$REML) "REML) = " else "ML) = ",
  #            format(x$logLik, digits=digits), "\n"))
  # the following adapted from print.summary.lme()
  dims <- x$dims
  cat("\nNumber of Observations:", dims[["N"]])
  cat("\nNumber of Groups: ")
  Ngrps <- dims$ngrps[1:dims$Q]
  if ((lNgrps <- length(Ngrps)) == 1) {
    cat(Ngrps, "\n")
  }
  else {
    sNgrps <- 1:lNgrps
    aux <- rep(names(Ngrps), sNgrps)
    aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps,
                                                   lNgrps))[!lower.tri(diag(lNgrps))])
    names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
    cat("\n")
    print(rev(Ngrps), ...)
  }
  print(x$fitstats)
  invisible(x)
}

S.glmerMod <- function(object, brief=FALSE, correlation=FALSE, exponentiate, ...){
  if (missing(exponentiate)) exponentiate <- object@resp$family$link %in% c("log", "logit")
  sumry <- summary(object)
  coefs <- sumry$coefficients
  # the following code for no of groups and obs borrowed from print.merMod()
  dims <- object@devcomp$dims
  ngrps <- vapply(object@flist, nlevels, 0L)
  vcov <- as.matrix(sumry$vcov)
  exp <- if (exponentiate) Confint(object, exponentiate=TRUE, silent=TRUE) else NULL
  result <- list(logLik=sumry$logLik, fixed.effects=coefs, random.effects=sumry$varcor,
                 call=sumry$call, brief=brief,
                 vcov=vcov, correlation=correlation, nobs=dims[["n"]], ngrps=ngrps,
                 exponentiate=exp, fitstats=fitstats(object))
  class(result) <- "S.glmerMod"
  result
}

print.S.glmerMod <- function(x, digits=max(3, getOption("digits") - 3),
                                  signif.stars = getOption("show.signif.stars"), ...){
  if (!x$brief) {
      cat("Generalized linear mixed model fit by ML\n")
      cat(formatCall(x$call))
  }
  cat("\nEstimates of Fixed Effects:\n")
  printCoefmat(x$fixed.effects, digits=digits, signif.stars=signif.stars)
  if (x$correlation) {
    # the following code adapted from print.summary.merMod()
    cor <- cov2cor(x$vcov)
    p <- ncol(cor)
    if (p > 1) {
      rn <- rownames(x$fixed.effects)
      rns <- abbreviate(rn, minlength = 11)
      cat("\nCorrelations of Fixed Effects:\n")
      cor <- matrix(format(round(cor, 3), nsmall = 3),
                    ncol = p, dimnames = list(rns, abbreviate(rn, minlength = 6)))
      cor[!lower.tri(cor)] <- ""
      print(cor[-1, -p, drop = FALSE], quote = FALSE)
    }
  }
  if (!is.null(x$exponentiate)){
    cat("\nExponentiated Fixed Effects and Confidence Bounds:\n")
    print(x$exponentiate)
  }
  cat("\nEstimates of Random Effects (Covariance Components):\n")
  print(x$random.effects, digits=digits)
  # cat(paste0("\nLog-likelihood = ", format(x$logLik, digits=digits), "\n"))
  # the following code adapted from lme4:::.prt.grps()
  cat(sprintf("\nNumber of obs: %d, groups: ", x$nobs),
      paste(paste(names(x$ngrps), x$ngrps, sep = ", "), collapse = "; "), fill = TRUE)
  print(x$fitstats)
  invisible(x)
}

Confint <- function(object, ...){
    UseMethod("Confint")
}

Confint.default <- function(object, estimate=TRUE, level=0.95, vcov., ...) {
    if (missing(vcov.)) result <- confint(object, level=level, ...)
    else{
        vc <- if (is.function(vcov.)) vcov.(object) else vcov.
        b <- coef(object)
        se <- sqrt(diag(vc))
        p <- 1 - (1 - level)/2
        z <- qnorm(p)
        result <- cbind(b - z*se, b + z*se)
        colnames(result) <- format.perc(c(1 - p, p), 3)
    }
    if (estimate){
        result <- cbind(coef(object), result)
        colnames(result)[1] <- "Estimate"
    }
    result
}

Confint.lm <- function(object, estimate=TRUE, parm, level = 0.95, vcov.= vcov(object, complete=FALSE), ...) {
    if (!missing(vcov.)) cat("Standard errors computed by", deparse(substitute(vcov.)), "\n")
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm))
        parm <- pnames
    else if (is.numeric(parm))
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    fac <- qt(a, object$df.residual)
    pct <- format.perc(a, 3)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm,  pct))
    ses <- sqrt(diag(if(is.matrix(vcov.)) vcov. else vcov.(object)))[parm]
    ci[] <- cf[parm] + ses %o% fac
    ci
    if (estimate){
        ci <- cbind(coef(object), ci)
        colnames(ci)[1] <- "Estimate"
    }
    ci
}

Confint.glm <- function(object, estimate=TRUE, exponentiate=FALSE, vcov., dispersion, type=c("LR", "Wald"), ...){
    type <- match.arg(type)
    silent <- list(...)$silent
    if (!missing(vcov.) && !missing(dispersion))
        stop("cannot specify both vcov. and dispersion arguments")
    if (!missing(vcov.) && (is.null(silent) || !silent)) cat("Standard errors computed by", deparse(substitute(vcov.)), "\n")
    result <- if (!missing(vcov.)) Confint.default(object, estimate=FALSE, vcov.=vcov(object, complete=FALSE), ...)
    else if (!missing(dispersion))
        Confint.default(object, estimate=FALSE, vcov.=dispersion*summary(object)$cov.unscaled, ...)
    else if (type == "LR"){
      suppressMessages(confint(object, ...))
    }
    else Confint.default(object, estimate=FALSE, vcov.=vcov(object))
    if (estimate){
        result <- cbind(coef(object), result)
        colnames(result)[1] <- "Estimate"
    }
    if (exponentiate){
        if (!object$family$link %in% c("log", "logit"))
            stop("exponentiated coefficients available only for log or logit link")
        if (is.null(silent) || !silent) cat("\nExponentiated Coefficients and Confidence Bounds\n")
        return(exp(result))
    }
    else return(result)
}

Confint.polr <- function(object, estimate=TRUE, exponentiate=FALSE, thresholds=!exponentiate, ...){
    dots <- list(...)
    level <- if (is.null(dots$level)) 0.95 else dots$level
    result <- suppressMessages(confint(object, ...))
    if (!is.matrix(result)) {
      cnames <- names(result)
      result <- matrix(result, nrow=1)
      colnames(result) <- cnames
      rownames(result) <- names(coef(object))
    }
    cnames <- colnames(result)
    if (estimate){
        result <- cbind(coef(object), result)
        colnames(result)[1] <- "Estimate"
    }
    if (thresholds){
        z <- qnorm(1 - (1 - level)/2)
        sumry <- suppressMessages(summary(object)$coefficients)
        sumry <- sumry[-(1:nrow(result)), ]
        b <- sumry[, 1]
        se <- sumry[, 2]
        sumry <- cbind(b - z*se, b + z*se)
        colnames(sumry) <- cnames
        if (estimate) {
            sumry <- cbind(b, sumry)
        }
        result <- rbind(result, sumry)
    }
    if (exponentiate) exp(result) else result
}

Confint.multinom <- function(object, estimate=TRUE, exponentiate=FALSE, ...){
    result <- confint(object)
    levs <- object$lev
    n.levs <- length(levs)
    b.names <- object$vcoefnames
    if (n.levs == 2){
        b <- coef(object)
        result <- cbind(b, result)
        colnames(result)[1] <- "Estimate"
        rownames(result) <- b.names
    }
    else if (estimate) {
        b <- object$wts
        b <- matrix(b, ncol=n.levs)
        b <- b[-1, , drop=FALSE]
        b <- b[ , -1, drop=FALSE]
        rownames(b) <- b.names
        colnames(b) <- levs[-1]
        result <- abind(b, result, along=2)
        dimnames(result)[[2]][1] <- "Estimate"
    }
    if (exponentiate) exp(result) else result
}

Confint.lme <- function(object, estimate=TRUE, level = 0.95, ...) {
  cf <- object$coefficients$fixed
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(cf), 2L), dimnames = list(names(cf),  pct))
  ses <- sqrt(diag(vcov(object, complete=FALSE)))
  ci[] <- cf + ses %o% fac
  if (estimate){
    ci <- cbind(cf, ci)
    colnames(ci)[1] <- "Estimate"
  }
  ci
}

Confint.lmerMod <- function(object, estimate=TRUE, level = 0.95, ...) {
  cf <- lme4::fixef(object)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(cf), 2L), dimnames = list(names(cf),  pct))
  ses <- sqrt(diag(as.matrix(vcov(object, complete=FALSE))))
  ci[] <- cf + ses %o% fac
  if (estimate){
    ci <- cbind(cf, ci)
    colnames(ci)[1] <- "Estimate"
  }
  ci
}

Confint.glmerMod <- function(object, estimate=TRUE, level = 0.95, exponentiate=FALSE, ...) {
  silent <- list(...)$silent
  cf <- lme4::fixef(object)
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qnorm(a)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(cf), 2L), dimnames = list(names(cf),  pct))
  ses <- sqrt(diag(as.matrix(vcov(object, complete=FALSE))))
  ci[] <- cf + ses %o% fac
  if (estimate){
    ci <- cbind(cf, ci)
    colnames(ci)[1] <- "Estimate"
  }
  if (exponentiate){
    if (!object@resp$family$link %in% c("log", "logit"))
      stop("exponentiated coefficients available only for log or logit link")
    if (is.null(silent) || !silent) cat("\nExponentiated Coefficients and Confidence Bounds\n")
    return(exp(ci))
  }
  else return(ci)
}

# the following function is not exported

coefs2use <- function(model, terms, intercept){
    vform <- update(formula(model), terms)
    if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
        stop("Only predictors in the formula can be used.")
    terms.model <- attr(attr(model.frame(model), "terms"), "term.labels")
    terms.vform <- attr(terms(vform), "term.labels")
    terms.used <- match(terms.vform, terms.model)
    mm <- model.matrix(model)
    model.names <- attributes(mm)$dimnames[[2]]
    model.assign <- attributes(mm)$assign
    use <- model.names[!is.na(match(model.assign, terms.used))]
    if (intercept && has.intercept(model)) c("(Intercept)", use) else use
}

# S <- function(model, terms, intercept, pvalues, digits, horizontal, ...){
#     UseMethod("S")
# }
#
#
# S.default <- function(model, terms = ~ ., intercept=missing(terms), pvalues=FALSE, digits=3, horizontal=TRUE, ...){
#     use <- coefs2use(model, terms, intercept)
#     sumry <- summary(model)
#     cols <- if (pvalues) c(1, 2, 4) else 1:2
#     coefs <- sumry$coefficients[use, cols, drop=FALSE]
#     colnames(coefs) <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|z|)") else c("Estimate", "Std. Error")
#     print(if (horizontal) t(coefs) else coefs, digits=digits)
#     invisible(sumry)
# }
#
# S.lm <- function(model, terms = ~ ., intercept=missing(terms), pvalues=FALSE, digits=3, horizontal=TRUE, vcov.=vcov, ...){
#     use <- coefs2use(model, terms, intercept)
#     sumry <- S(model, vcov.=vcov., ...)
#     cols <- if (pvalues) c(1, 2, 4) else 1:2
#     coefs <- sumry$coefficients[use, cols, drop=FALSE]
#     colnames(coefs) <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|t|)") else c("Estimate", "Std. Error")
#     print(if (horizontal) t(coefs) else coefs, digits=digits)
#     if (missing(terms)) cat("\n Residual SD =", format(sumry$sigma, digits=digits),
#                             "on", model$df.residual, "df, R-squared =", format(sumry$r.squared, digits=digits))
#     invisible(sumry)
# }
#
# S.glm <- function(model, terms = ~ ., intercept=missing(terms), pvalues=FALSE, digits=3, horizontal=TRUE, vcov., dispersion, exponentiate, ...){
#     if (!missing(vcov.) && !missing(dispersion))
#         stop("cannot specify both the dispersion and vcov. arguments")
#     if (missing(exponentiate)) exponentiate <- model$family$link %in% c("log", "logit")
#     use <- coefs2use(model, terms, intercept)
#     sumry <- if (!missing(vcov.)) S(model, digits, vcov.=vcov., ...)
#     else if (!missing(dispersion)) S(model, digits, dispersion=dispersion, ...)
#     else summary(model, ...)
#     cols <- if (pvalues) c(1, 2, 4) else 1:2
#     coefs <- sumry$coefficients[use, cols, drop=FALSE]
#     colnames(coefs) <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|z|)") else c("Estimate", "Std. Error")
#     if (exponentiate){
#         coefs <- cbind(coefs, exp(coefs[, 1]))
#         colnames(coefs)[if (pvalues) 4 else 3] <- "exp(Estimate)"
#     }
#     print(if (horizontal) t(coefs) else coefs, digits=digits)
#     if (missing(terms)) cat("\n Residual deviance =", format(model$deviance, digits=digits),
#                             "on", model$df.residual, "df",
#                             if (family(model)$family %in% c("binomial", "poisson")) ""
#                             else (paste(", Est. dispersion =", format(sumry$dispersion, digits=digits))))
#     invisible(sumry)
# }
#
# S.polr <- function(model, terms = ~ ., intercept, pvalues=FALSE, digits=3, horizontal=TRUE, exponentiate=TRUE, ...){
#     sumry <- summary(model)
#     coefs <- sumry$coefficients[ , 1:2]
#     if (pvalues) {
#         coefs <- cbind(coefs, 2*pnorm(abs(coefs[ , 1]/coefs[, 2]), lower.tail=FALSE))
#     }
#     use <- if (missing(terms)) 1:nrow(coefs) else coefs2use(model, terms, FALSE)
#     coefs <- coefs[use, , drop=FALSE]
#     colnames(coefs) <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|z|)") else c("Estimate", "Std. Error")
#     if (exponentiate){
#         coefs <- cbind(coefs, exp(coefs[, 1]))
#         colnames(coefs)[if (pvalues) 4 else 3] <- "exp(Estimate)"
#         if (missing(terms)){
#             n.thresholds <- length(model$zeta)
#             n.pars <- nrow(coefs)
#             coefs[(n.pars - n.thresholds + 1):n.pars , if (pvalues) 4 else 3] <- NA
#         }
#     }
#     print(if (horizontal) t(coefs) else coefs, digits=digits, na.print="")
#     if (missing(terms)) cat("\n Residual deviance =", format(model$deviance, digits=digits),
#         "on", model$df.residual, "df")
#     invisible(sumry)
# }
#
# S.multinom <- function(model, terms = ~ ., intercept=missing(terms), pvalues=FALSE, digits=3, horizontal=TRUE, exponentiate=TRUE, ...){
#     use <- coefs2use(model, terms, intercept)
#     sumry <- summary(model, ...)
#     b <- sumry$coefficients
#     se <- sumry$standard.errors
#     p <- 2*pnorm(abs(b/se), lower.tail=FALSE)
#     levels <- sumry$lev
#     labels <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|z|)") else c("Estimate", "Std. Error")
#     if (exponentiate) labels <- c(labels, "exp(Estimate)")
#     if (length(levels) == 2){
#         b <- b[use]
#         se <- se[use]
#         p <- p[use]
#         table <- if (pvalues) rbind(b, se, p) else rbind(b, se)
#         if (exponentiate) table <- rbind(table, exp(b))
#         rownames(table) <- labels
#         cat("\n ", levels[2], "\n")
#         print(if (horizontal) table else t(table), digits=digits)
#     }
#     else{
#         b <- b[, use, drop=FALSE]
#         se <- se[, use, drop=FALSE]
#         p <- p[, use, drop=FALSE]
#         table <- if (pvalues) abind(t(b), t(se), t(p), along=1.5) else abind(t(b), t(se), along=1.5)
#         if (exponentiate) table <- abind(table, t(exp(b)), along=2)
#         dimnames(table)[[2]] <- labels
#         for (level in levels[-1]){
#             cat("\n ", level, "\n")
#             result <- if (horizontal) t(table[, , level]) else table[, , level]
#             if (dim(table)[1] == 1){
#                 if (horizontal) rownames(result) <- dimnames(table)[1] else {
#                     result <- matrix(result, ncol=1)
#                     colnames(result) <- dimnames(table)[1]
#                 }
#             }
#             print(result, digits=digits)
#         }
#     }
#     if (missing(terms)) cat("\n Residual deviance =", format(model$deviance, digits=digits),
#         "fitting", length(b), "parameters")
#     invisible(sumry)
# }

