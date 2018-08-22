# added by J. Fox on 2017-10-14

poTest <- function(model, ...){
    UseMethod("poTest")
}

poTest.polr <- function(model, ...){
    if (model$method != "logistic") stop("test for proportional odds is only for the logistic model")
    X <- model.matrix(model)
    y <- model.frame(model)[, 1]
    levels <- levels(y)
    k <- length(levels)
    p <- ncol(X) - 1
    y <- as.numeric(y)
    models <- vector(k - 1, mode="list")
    for (j in 1:(k - 1)){
        models[[j]] <- glm(y > j ~ X - 1, family=binomial)
    }
    vcov <- matrix(0, (k - 1)*p, (k - 1)*p)
    for (el in 1:(k - 1)){
        for (j in 1:el){
            W.j.el <- fitted(models[[el]]) - fitted(models[[j]])*fitted(models[[el]])
            W.el.el <- fitted(models[[el]]) - fitted(models[[el]])^2
            W.j.j <- fitted(models[[j]]) - fitted(models[[j]])^2
            V <- solve(t(X * W.j.j) %*% X) %*% (t(X * W.j.el) %*% X) %*% solve(t(X * W.el.el) %*% X)
            subs.j <- (j - 1)*p + 1:p
            subs.el <- (el - 1)*p + 1:p
            vcov[subs.j, subs.el] <- vcov[subs.el, subs.j] <- V[-1, -1]
        }
    }
    beta <- unlist(lapply(models, function(m) coef(m)[-1]))
    D <- matrix(0, (k - 2)*p, (k - 1)*p)
    I <- diag(p)
    for (j in 1:(k - 2)){
        subs.j <- (j - 1)*p + 1:p
        subs.el <- j*p + 1:p
        D[subs.j, 1:p] <- I
        D[subs.j, subs.el] <- -I
    }
    chisq <- t(D %*% beta) %*% solve(D %*% vcov %*% t(D)) %*% (D %*% beta)
    df <- (k - 2)*p
    chisq.p <- numeric(p)
    zeros <- matrix(0, k - 2, (k - 1)*p)
    D.p <- vector(p, mode="list")
    for (i in 1:p){
        DD <- zeros
        j <- 1:(k - 2)
        DD[j, i] <- 1
        DD[cbind(j, j*p + i)] <- -1
        chisq.p[i] <- t(DD %*% beta) %*% solve(DD %*% vcov %*% t(DD)) %*% (DD %*% beta)
        D.p[[i]] <- DD
    }
    b <- coef(model)
    coef.names <- names(b)
    b <- cbind(b, matrix(beta, ncol = k - 1))
    colnames(b) <- c("b[polr]", paste0("b[>", levels[1:(k - 1)], "]"))
    result <- list(call=model$call, coef.names=coef.names, b=b,
                   vcov=vcov, D=D, chisq=as.vector(chisq), df=df,
                   D.p=D.p, chisq.p=chisq.p, df.p = k - 2)
    class(result) <- "poTest"
    result
}

print.poTest <- function(x, digits=3, ...){
    cat("\nTests for Proportional Odds\n")
    print(x$call)
    cat("\n")
    names <- c("Overall", x$coef.names)
    chisq <- c(x$chisq, x$chisq.p)
    df <- c(x$df, rep(x$df.p, length(x$chisq.p)))
    pval <- pchisq(chisq, df, lower.tail=FALSE)
    table <- cbind(chisq, df, pval)
    colnames(table) <- c("Chisquare", "df", "Pr(>Chisq)")
    b <- x$b
    b <- rbind(rep(NA, ncol(b)), b)
    table <- cbind(b, table)
    rownames(table) <- names
    printCoefmat(table, P.values=TRUE, has.Pvalue=TRUE, tst.ind = ncol(b) + 1,
                 na.print="", digits=digits)
    invisible(x)
}
