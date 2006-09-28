# last modified 28 Sept 2006 by J. Fox

linear.hypothesis <- function (model, ...)       
    UseMethod("linear.hypothesis")
    
lht <- function (model, ...)
    UseMethod("linear.hypothesis")

linear.hypothesis.default <- function(model, hypothesis.matrix, rhs=NULL, 
    test=c("Chisq", "F"), vcov.=NULL, verbose=FALSE, ...){
    makeHypothesis <- function(cnames, hypothesis, rhs = NULL){
        parseTerms <- function(terms){
            component <- gsub("^[-\\ 0-9\\.]+", "", terms)
            component <- gsub(" ", "", component, extended=FALSE, fixed=TRUE)
            component
            }
        stripchars <- function(x) {
            x <- gsub(" ", "", x, extended = FALSE, fixed = TRUE)
            x <- gsub("*", "", x, extended = FALSE, fixed = TRUE)
            x <- gsub("-", "+-", x, extended = FALSE, fixed = TRUE)
            x <- strsplit(x, "+", extended = FALSE, fixed = TRUE)[[1]]
            x <- x[x!=""]
            x
            }
        char2num <- function(x) {
            x[x == ""] <- "1"
            x[x == "-"] <- "-1"
            as.numeric(x)
            }
        constants <- function(x) {
            with.coef <- unique(unlist(sapply(cnames, 
                function(y) which(y == parseTerms(x)))))
            if(length(with.coef) > 0) x <- x[-with.coef]
            x <- if(is.null(x)) 0 else sum(as.numeric(x))
            if(any(is.na(x)))
              stop('The hypothesis "', hypothesis, 
                '" is not well formed: contains bad coefficient/variable names.')
            x
            }
        coefvector <- function(x, y) {
            rv <- gsub(" ", "", x, extended=FALSE, fixed=TRUE) ==
                parseTerms(y)
            if (!any(rv)) return(0)
            if (sum(rv) > 1) stop('The hypothesis "', hypothesis, 
                '" is not well formed.')
            rv <- sum(char2num(unlist(strsplit(y[rv], x, extended=FALSE, fixed=TRUE))))
            if(is.na(rv))
              stop('The hypothesis "', hypothesis, 
                '" is not well formed: contains non-numeric coefficients.')
            rv
            }
        rhs <- rep(rhs, length.out = length(hypothesis))
        if(length(hypothesis) > 1)
            return(rbind(Recall(cnames, hypothesis[1], rhs[1]), 
                Recall(cnames, hypothesis[-1], rhs[-1])))
        lhs <- strsplit(hypothesis, "=", extended = FALSE, fixed = TRUE)[[1]]
        if(is.null(rhs)) {
            if(length(lhs) < 2) rhs <- "0"
                else if(length(lhs) == 2) {
                    rhs <- lhs[2]
                    lhs <- lhs[1]
                    } 
                else stop('The hypothesis "', hypothesis, 
                    '" is not well formed: contains more than one = sign.')
                } 
            else {
                if(length(lhs) < 2) as.character(rhs)
                else stop('The hypothesis "', hypothesis,
                    '" is not well formed: contains a = sign although rhs was specified.')
                }
        lhs <- stripchars(lhs)
        rhs <- stripchars(rhs)
        rval <- sapply(cnames, coefvector, y = lhs) - sapply(cnames, coefvector, y = rhs)
        rval <- c(rval, constants(rhs) - constants(lhs))
        names(rval) <- c(cnames, "*rhs*")  
        rval
        }
    printHypothesis <- function(L, rhs, cnames){
        hyp <- rep("", nrow(L))
        for (i in 1:nrow(L)){
            sel <- L[i,] != 0
            h <- L[i, sel]
            h <- ifelse(h < 0, as.character(h), paste("+", h, sep="")) 
            nms <- cnames[sel]
            h <- paste(h, nms) 
            h <- gsub("-1", "-", h)
            h <- gsub("+1", "+", h, extended=FALSE, fixed=TRUE)
            h <- gsub("-", " - ", h)
            h <- gsub("+", "  + ", h, extended=FALSE, fixed=TRUE)
            h <- paste(h, collapse="")
            h <- gsub("  ", " ", h, extended=FALSE, fixed=TRUE)
            h <- sub("^\\ \\+", "", h)
            h <- sub("^\\ ", "", h)
            h <- sub("^-\\ ", "-", h)
            hyp[i] <- paste(h, "=", rhs[i])
            }
        hyp
        }
    df <- df.residual(model)
    if(is.null(df)) df <- Inf ## if no residual df available
    V <- if(is.null(vcov.)) vcov(model)  
        else if(is.function(vcov.)) vcov.(model) else vcov.
    b <- coef(model)
    if(is.character(hypothesis.matrix)) {    
        L <- makeHypothesis(names(b), hypothesis.matrix, rhs)
        if(is.null(dim(L))) L <- t(L)
        rhs <- L[, NCOL(L)]
        L <- L[, -NCOL(L), drop = FALSE]
        rownames(L) <- hypothesis.matrix
        } 
    else {
        L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix) 
        else hypothesis.matrix
        if(is.null(rhs)) rhs <- rep(0, nrow(L))
        }  
    q <- NROW(L)
    if (verbose){
        cat("\nHypothesis matrix:\n")    
        print(L)
        cat("\nRight-hand-side vector:\n")
        print(rhs)
        cat("\n")
        }
    SSH <- as.vector(t(L %*% b - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% b - rhs))
    test <- match.arg(test)
    if(!(is.finite(df) && df > 0)) test <- "Chisq"
    name <- try(formula(model), silent = TRUE)
    if(inherits(name, "try-error")) name <- substitute(model)  
    title <- "Linear hypothesis test\n\nHypothesis:"
    topnote <- paste("Model 1: ", paste(deparse(name), collapse = "\n"), "\n",
    	   "Model 2: restricted model", sep = "")
    note <- if (is.null(vcov.)) "" 
        else "\nNote: Coefficient covariance matrix supplied.\n"
    rval <- matrix(rep(NA, 8), ncol = 4)
    colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
    rownames(rval) <- 1:2
    rval[,1] <- c(df, df+q) 
    if(test == "F") {
        f <- SSH/q
        p <- pf(f, q, df, lower.tail = FALSE)
        rval[2,2:4] <- c(-q, f, p)
        } 
    else {
        p <- pchisq(SSH, q, lower.tail = FALSE)
        rval[2,2:4] <- c(-q, SSH, p)
        }
    if(!(is.finite(df) && df > 0)) rval <- rval[,-1]
    structure(as.data.frame(rval), 
        heading = c(title, printHypothesis(L, rhs, names(b)), "", topnote, note), 
        class = c("anova", "data.frame"))
    }

linear.hypothesis.glm <- function(model, ...)
    linear.hypothesis.default(model, ...)

linear.hypothesis.lm <- function(model, hypothesis.matrix, rhs=NULL,
    test=c("F", "Chisq"), vcov.=NULL, white.adjust=FALSE, ...){
    if (is.aliased(model)) stop("One or more terms aliased in model.")
    test <- match.arg(test)
    if(identical(white.adjust, TRUE)) white.adjust <- "hc3"
    if(is.null(vcov.) && is.character(white.adjust))
        vcov. <- hccm(model, type = white.adjust)
    rval <- linear.hypothesis.default(model, hypothesis.matrix, rhs = rhs,
        test = test, vcov. = vcov., ...)
    if(is.null(vcov.)) {
        rval2 <- matrix(rep(NA, 4), ncol = 2)
        colnames(rval2) <- c("RSS", "Sum of Sq")
        SSH <- rval[2,test]
        if(test == "F") SSH <- SSH * abs(rval[2, "Df"])
        df <- rval[1, "Res.Df"]
        error.SS <- deviance(model)
        rval2[,1] <- c(error.SS, error.SS + SSH * error.SS/df)
        rval2[2,2] <- -diff(rval2[,1])
        rval2 <- cbind(rval, rval2)[,c(1, 5, 2, 6, 3, 4)]
        class(rval2) <- c("anova", "data.frame")
        attr(rval2, "heading") <- attr(rval, "heading")
        rval <- rval2
        }
    rval
    }



