# last modified 27 Dec 2008 by J. Fox

vcov.default <- function(object, ...){
	stop(paste("there is no vcov() method for models of class", 
			paste(class(object), collapse=", ")))
}

has.intercept.matrix <- function (model, ...) {
    "(Intercept)" %in% colnames(model)
    }
    

makeHypothesis <- function(cnames, hypothesis, rhs = NULL){
    parseTerms <- function(terms){
        component <- gsub("^[-\\ 0-9\\.]+", "", terms)
        component <- gsub(" ", "", component, fixed=TRUE)
        component
        }
    stripchars <- function(x) {
        x <- gsub(" ", "", x, fixed = TRUE)
        x <- gsub("*", "", x, fixed = TRUE)
        x <- gsub("-", "+-", x, fixed = TRUE)
        x <- strsplit(x, "+", fixed = TRUE)[[1]]
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
        if (length(with.coef) > 0) x <- x[-with.coef]
        x <- if (is.null(x)) 0 else sum(as.numeric(x))
        if (any(is.na(x)))
          stop('The hypothesis "', hypothesis, 
            '" is not well formed: contains bad coefficient/variable names.')
        x
        }
    coefvector <- function(x, y) {
        rv <- gsub(" ", "", x, fixed=TRUE) ==
            parseTerms(y)
        if (!any(rv)) return(0)
        if (sum(rv) > 1) stop('The hypothesis "', hypothesis, 
            '" is not well formed.')
        rv <- sum(char2num(unlist(strsplit(y[rv], x, fixed=TRUE))))
        if (is.na(rv))
          stop('The hypothesis "', hypothesis, 
            '" is not well formed: contains non-numeric coefficients.')
        rv
        }
    rhs <- rep(rhs, length.out = length(hypothesis))
    if (length(hypothesis) > 1)
        return(rbind(Recall(cnames, hypothesis[1], rhs[1]), 
            Recall(cnames, hypothesis[-1], rhs[-1])))
    lhs <- strsplit(hypothesis, "=", fixed=TRUE)[[1]]
    if (is.null(rhs)) {
        if (length(lhs) < 2) rhs <- "0"
            else if (length(lhs) == 2) {
                rhs <- lhs[2]
                lhs <- lhs[1]
                } 
            else stop('The hypothesis "', hypothesis, 
                '" is not well formed: contains more than one = sign.')
            } 
        else {
            if (length(lhs) < 2) as.character(rhs)
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
        h <- gsub("+1", "+", h, fixed=TRUE)
        h <- gsub("-", " - ", h)
        h <- gsub("+", "  + ", h, fixed=TRUE)
        h <- paste(h, collapse="")
        h <- gsub("  ", " ", h, fixed=TRUE)
        h <- sub("^\\ \\+", "", h)
        h <- sub("^\\ ", "", h)
        h <- sub("^-\\ ", "-", h)
        hyp[i] <- paste(h, "=", rhs[i])
        }
    hyp
    }
    
linear.hypothesis <- function (model, ...)       
    UseMethod("linear.hypothesis")
    
lht <- function (model, ...)
    UseMethod("linear.hypothesis")

linear.hypothesis.default <- function(model, hypothesis.matrix, rhs=NULL, 
    test=c("Chisq", "F"), vcov.=NULL, verbose=FALSE, ...){
    df <- df.residual(model)
    if (is.null(df)) df <- Inf ## if no residual df available
    V <- if (is.null(vcov.)) vcov(model)  
        else if (is.function(vcov.)) vcov.(model) else vcov.
    b <- coef(model)
	if (is.null(b)) stop(paste("there is no coef() method for models of class", 
				paste(class(model), collapse=", ")))
    if (is.character(hypothesis.matrix)) {    
        L <- makeHypothesis(names(b), hypothesis.matrix, rhs)
        if (is.null(dim(L))) L <- t(L)
        rhs <- L[, NCOL(L)]
        L <- L[, -NCOL(L), drop = FALSE]
        rownames(L) <- hypothesis.matrix
        } 
    else {
        L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix) 
        else hypothesis.matrix
        if (is.null(rhs)) rhs <- rep(0, nrow(L))
        }  
    q <- NROW(L)
    if (verbose){
        cat("\nHypothesis matrix:\n")    
        print(L)
        cat("\nRight-hand-side vector:\n")
        print(rhs)
		cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs)\n")
		print(drop(L %*% b - rhs))
        cat("\n")
        }
    SSH <- as.vector(t(L %*% b - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% b - rhs))
    test <- match.arg(test)
    if (!(is.finite(df) && df > 0)) test <- "Chisq"
    name <- try(formula(model), silent = TRUE)
    if (inherits(name, "try-error")) name <- substitute(model)  
    title <- "Linear hypothesis test\n\nHypothesis:"
    topnote <- paste("Model 1: ", paste(deparse(name), collapse = "\n"), "\n",
    	   "Model 2: restricted model", sep = "")
    note <- if (is.null(vcov.)) "" 
        else "\nNote: Coefficient covariance matrix supplied.\n"
    rval <- matrix(rep(NA, 8), ncol = 4)
    colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
    rownames(rval) <- 1:2
    rval[,1] <- c(df, df+q) 
    if (test == "F") {
        f <- SSH/q
        p <- pf(f, q, df, lower.tail = FALSE)
        rval[2,2:4] <- c(-q, f, p)
        } 
    else {
        p <- pchisq(SSH, q, lower.tail = FALSE)
        rval[2,2:4] <- c(-q, SSH, p)
        }
    if (!(is.finite(df) && df > 0)) rval <- rval[,-1]
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
    if (identical(white.adjust, TRUE)) white.adjust <- "hc3"
    if (is.null(vcov.) && is.character(white.adjust))
        vcov. <- hccm(model, type = white.adjust)
    rval <- linear.hypothesis.default(model, hypothesis.matrix, rhs = rhs,
        test = test, vcov. = vcov., ...)
    if (is.null(vcov.)) {
        rval2 <- matrix(rep(NA, 4), ncol = 2)
        colnames(rval2) <- c("RSS", "Sum of Sq")
        SSH <- rval[2,test]
        if (test == "F") SSH <- SSH * abs(rval[2, "Df"])
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
    
linear.hypothesis.mlm <- function(model, hypothesis.matrix, rhs=NULL, SSPE, V,
	test, idata, icontrasts=c("contr.sum", "contr.poly"), idesign, iterms, 
	P=NULL, title="", verbose=FALSE, ...){
	check <- function(X){ # check block orthogonality of model matrix
		XX <- crossprod(X)
		terms <- attr(X, "assign")
		for (term in unique(terms)){
			subs <- term == terms
			XX[subs, subs] <- 0
		}
		!any(abs(XX) > sqrt(.Machine$double.eps))
	}
	if (missing(test)) test <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
	test <- match.arg(test, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
		several.ok=TRUE)              
	df.residual <- df.residual(model)
	if (missing (V)) V <- solve(crossprod(model.matrix(model)))
	B <- coef(model)
	if (is.character(hypothesis.matrix)) {    
		L <- makeHypothesis(rownames(B), hypothesis.matrix, rhs)
		if (is.null(dim(L))) L <- t(L)
		L <- L[, -NCOL(L), drop = FALSE]
		rownames(L) <- hypothesis.matrix
	} 
	else {
		L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix) 
			else hypothesis.matrix
	}
	if (missing(SSPE)) SSPE <- crossprod(residuals(model))
	if (!missing(idata)){
		for (i in 1:length(idata)){
			if (is.null(attr(idata[,i], "contrasts"))){
				contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
					else icontrasts[1]
			}
		}
		if (missing(idesign)) stop("idesign (intra-subject design) missing.")
		X.design <- model.matrix(idesign, data=idata)
		if (!check(X.design)) 
			stop("Terms in the intra-subject model matrix are not orthogonal.")
		intercept <- has.intercept(X.design)
		term.names <- term.names(idesign)
		if (intercept) term.names <- c("(Intercept)", term.names)
		which.terms <- match(iterms, term.names)
		if (any(nas <- is.na(which.terms))){
			if (sum(nas) == 1) 
				stop('The term "', iterms[nas],'" is not in the intrasubject design.')
			else stop("The following terms are not in the intrasubject design: ",
					paste(iterms[nas], collapse=", "), ".")
		}
		select <- apply(outer(which.terms, attr(X.design, "assign") + intercept, "=="),
			2, any)
		P <- X.design[, select, drop=FALSE]
	}
	if (!is.null(P)){
		rownames(P) <- colnames(B)   
		SSPE <- t(P) %*% SSPE %*% P
		B <- B %*% P
	}
	rank <- sum(eigen(SSPE, only.values=TRUE)$values >= sqrt(.Machine$double.eps))
	if (rank < ncol(SSPE)) 
		stop("The error SSP matrix is apparently of deficient rank = ",
			rank, " < ", ncol(SSPE))
	r <- ncol(B)
	if (is.null(rhs)) rhs <- matrix(0, nrow(L), r)
	rownames(rhs) <- rownames(L)
	colnames(rhs) <- colnames(B)
	q <- NROW(L)
	if (verbose){
		cat("\nHypothesis matrix:\n")    
		print(L)
		cat("\nRight-hand-side matrix:\n")
		print(rhs)
		cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs):\n")
		print(drop(L %*% B - rhs))
		cat("\n")
	}
	SSPH <- t(L %*% B - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% B - rhs)
	rval <- list(SSPH=SSPH, SSPE=SSPE, df=q, r=r, df.residual=df.residual, P=P,
		title=title, test=test)
	class(rval) <- "linear.hypothesis.mlm"
	rval
}
    
print.linear.hypothesis.mlm <- function(x, SSP=TRUE, SSPE=SSP, 
    digits=unlist(options("digits")), ...){
    test <- x$test
    if (!is.null(x$P) && SSP){
        P <- x$P
        cat("\n Response transformation matrix:\n")
        attr(P, "assign") <- NULL
        attr(P, "contrasts") <- NULL
        print(P, digits=digits)
        } 
    if (SSP){
        cat("\nSum of squares and products for the hypothesis:\n")
        print(x$SSPH, digits=digits)
        }
    if (SSPE){
        cat("\nSum of squares and products for error:\n")
        print(x$SSPE, digits=digits)
        }     
    SSPE.qr <- qr(x$SSPE)    
    # the following code is adapted from summary.manova
    eigs <- Re(eigen(qr.coef(SSPE.qr, x$SSPH), symmetric = FALSE)$values)
    tests <- matrix(NA, 4, 4)
    rownames(tests) <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")  
    if ("Pillai" %in% test) 
        tests[1, 1:4] <- stats:::Pillai(eigs, x$df, x$df.residual)
    if ("Wilks" %in% test) 
        tests[2, 1:4] <- stats:::Wilks(eigs, x$df, x$df.residual)
    if ("Hotelling-Lawley" %in% test)
        tests[3, 1:4] <- stats:::HL(eigs, x$df, x$df.residual)
    if ("Roy" %in% test)
        tests[4, 1:4] <- stats:::Roy(eigs, x$df, x$df.residual)
    tests <- na.omit(tests)
    ok <- tests[, 2] >= 0 & tests[, 3] > 0 & tests[, 4] > 0
    ok <- !is.na(ok) & ok
    tests <- cbind(x$df, tests, pf(tests[ok, 2], tests[ok, 3], tests[ok, 4], 
        lower.tail = FALSE))
    colnames(tests) <- c("Df", "test stat", "approx F", "num Df", "den Df", "Pr(>F)")  
    tests <- structure(as.data.frame(tests), 
        heading = paste("\nMultivariate Test",
            if (nrow(tests) > 1) "s", ": ", x$title, sep=""), 
        class = c("anova", "data.frame"))
    print(tests, digits=digits)
    invisible(x)
    }

linear.hypothesis.survreg <- function(model, hypothesis.matrix, rhs=NULL, 
	test=c("Chisq", "F"), vcov., verbose=FALSE, ...){
	if (missing(vcov.)) {
		vcov. <- vcov(model)
		p <- nrow(vcov.)
		vcov. <- vcov.[-p, -p]
	}
	linear.hypothesis.default(model, hypothesis.matrix, rhs, test, vcov.)
}