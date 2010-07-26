#-------------------------------------------------------------------------------
# Revision history:
#   2009-01-16: replaced unlist(options("foo")) with getOption("foo")
#   2009-09-16: optionally allow models with aliased coefficients. J. Fox
#   2009-12-10: modification by A. Zeileis to allow wider range of coef. names.
#   2009-12-22: small changes to linearHypothesis.mlm() to handle user-specified
#               within-subjects designs in Anova()
#   2010-05-21: linearHypothesis.default() and .lm() changed so that differences
#               in df, etc. will be postive.
#   2010-06-12: linearHypothesis.mlm() changed to allow observation weights
#	2010-06-22: fixed bug in linearHypothesis.lm caused by 2010-05-21 revision
#-------------------------------------------------------------------------------

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
	constants <- function(x, y) { 
		with.coef <- unique(unlist(sapply(y,
								function(z) which(z == parseTerms(x)))))
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
	
	cnames_symb <- sapply(c("@", "#", "~"), function(x) length(grep(x, cnames)) < 1)
	
	if(any(cnames_symb)) {
		cnames_symb <- head(c("@", "#", "~")[cnames_symb], 1)
		cnames_symb <- paste(cnames_symb, seq_along(cnames), cnames_symb, sep = "")
		hypothesis_symb <- hypothesis
		for(i in order(nchar(cnames), decreasing = TRUE))
			hypothesis_symb <- gsub(cnames[i], cnames_symb[i], hypothesis_symb, fixed = TRUE)
	} else {
		stop('The hypothesis "', hypothesis,
				'" is not well formed: contains non-standard coefficient names.')
	}
	
	lhs <- strsplit(hypothesis_symb, "=", fixed=TRUE)[[1]] 
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
	rval <- sapply(cnames_symb, coefvector, y = lhs) - sapply(cnames_symb, coefvector, y = rhs) 
	rval <- c(rval, constants(rhs, cnames_symb) - constants(lhs, cnames_symb)) 
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

linearHypothesis <- function (model, ...)
	UseMethod("linearHypothesis")

lht <- function (model, ...)
	UseMethod("linearHypothesis")

linearHypothesis.default <- function(model, hypothesis.matrix, rhs=NULL,
		test=c("Chisq", "F"), vcov.=NULL, singular.ok=FALSE, verbose=FALSE, ...){
	df <- df.residual(model)
	if (is.null(df)) df <- Inf ## if no residual df available
	V <- if (is.null(vcov.)) vcov(model)
			else if (is.function(vcov.)) vcov.(model) else vcov.
	b <- coef(model)
	if (any(aliased <- is.na(b)) && !singular.ok)
		stop("there are aliased coefficients in the model")
	b <- b[!aliased]
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
	topnote <- paste("Model 1: restricted model","\n", "Model 2: ", 
			paste(deparse(name), collapse = "\n"), sep = "")
	note <- if (is.null(vcov.)) ""
			else "\nNote: Coefficient covariance matrix supplied.\n"
	rval <- matrix(rep(NA, 8), ncol = 4)
	colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
	rownames(rval) <- 1:2
	rval[,1] <- c(df+q, df)
	if (test == "F") {
		f <- SSH/q
		p <- pf(f, q, df, lower.tail = FALSE)
		rval[2, 2:4] <- c(q, f, p)
	}
	else {
		p <- pchisq(SSH, q, lower.tail = FALSE)
		rval[2, 2:4] <- c(q, SSH, p)
	}
	if (!(is.finite(df) && df > 0)) rval <- rval[,-1]
	structure(as.data.frame(rval),
			heading = c(title, printHypothesis(L, rhs, names(b)), "", topnote, note),
			class = c("anova", "data.frame"))
}

linearHypothesis.glm <- function(model, ...)
	linearHypothesis.default(model, ...)

linearHypothesis.lm <- function(model, hypothesis.matrix, rhs=NULL,
		test=c("F", "Chisq"), vcov.=NULL,
		white.adjust=c(FALSE, TRUE, "hc3", "hc0", "hc1", "hc2", "hc4"),
		singular.ok=FALSE, ...){
	if (!singular.ok && is.aliased(model))
		stop("there are aliased coefficients in the model.")
	test <- match.arg(test)
	white.adjust <- as.character(white.adjust)
	white.adjust <- match.arg(white.adjust)
	if (white.adjust != "FALSE"){
		if (white.adjust == "TRUE") white.adjust <- "hc3"
		vcov. <- hccm(model, type=white.adjust)
	}
	rval <- linearHypothesis.default(model, hypothesis.matrix, rhs = rhs,
			test = test, vcov. = vcov., singular.ok=singular.ok, ...)
	if (is.null(vcov.)) {
		rval2 <- matrix(rep(NA, 4), ncol = 2)
		colnames(rval2) <- c("RSS", "Sum of Sq")
		SSH <- rval[2,test]
		if (test == "F") SSH <- SSH * abs(rval[2, "Df"])
		df <- rval[2, "Res.Df"]
		error.SS <- deviance(model)
		rval2[,1] <- c(error.SS + SSH * error.SS/df, error.SS)
		rval2[2,2] <- abs(diff(rval2[,1]))
		rval2 <- cbind(rval, rval2)[,c(1, 5, 2, 6, 3, 4)]
		class(rval2) <- c("anova", "data.frame")
		attr(rval2, "heading") <- attr(rval, "heading")
		rval <- rval2
	}
	rval
}


check.imatrix <- function(X, terms){ 
# check block orthogonality of within-subjects model matrix
	XX <- crossprod(X)
	if (missing(terms)) terms <- attr(X, "assign")
	for (term in unique(terms)){
		subs <- term == terms
		XX[subs, subs] <- 0
	}
	if (any(abs(XX) > sqrt(.Machine$double.eps)))
		stop("Terms in the intra-subject model matrix are not orthogonal.")
}

linearHypothesis.mlm <- function(model, hypothesis.matrix, rhs=NULL, SSPE, V,
		test, idata, icontrasts=c("contr.sum", "contr.poly"), idesign, iterms,
		check.imatrix=TRUE, P=NULL, title="", verbose=FALSE, ...){
	if (missing(test)) test <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
	test <- match.arg(test, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
			several.ok=TRUE)
	df.residual <- df.residual(model)
	wts <- if (!is.null(model$weights)) model$weights else rep(1,nrow(model.matrix(model)))
	# V = (X'WX)^{-1}
	if (missing (V)) V <- solve(wcrossprod(model.matrix(model), w=wts))
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
	# SSPE = E'WE
	if (missing(SSPE)) SSPE <- wcrossprod(residuals(model),w=wts)
	if (missing(idata)) idata <- NULL
	if (missing(idesign)) idesign <- NULL
	if (!is.null(idata)){
		for (i in 1:length(idata)){
			if (is.null(attr(idata[,i], "contrasts"))){
				contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
						else icontrasts[1]
			}
		}
		if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
		X.design <- model.matrix(idesign, data=idata)
		if (check.imatrix) check.imatrix(X.design)
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
	class(rval) <- "linearHypothesis.mlm"
	rval
}


#linearHypothesis.mlm <- function(model, hypothesis.matrix, rhs=NULL, SSPE, V,
#   test, idata, icontrasts=c("contr.sum", "contr.poly"), idesign, iterms,
#   check.imatrix=TRUE, P=NULL, title="", verbose=FALSE, ...){
#   if (missing(test)) test <- c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")
#   test <- match.arg(test, c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
#       several.ok=TRUE)
#   df.residual <- df.residual(model)
#   if (missing (V)) V <- solve(crossprod(model.matrix(model)))
#   B <- coef(model)
#   if (is.character(hypothesis.matrix)) {
#       L <- makeHypothesis(rownames(B), hypothesis.matrix, rhs)
#       if (is.null(dim(L))) L <- t(L)
#       L <- L[, -NCOL(L), drop = FALSE]
#       rownames(L) <- hypothesis.matrix
#   }
#   else {
#       L <- if (is.null(dim(hypothesis.matrix))) t(hypothesis.matrix)
#           else hypothesis.matrix
#   }
#   if (missing(SSPE)) SSPE <- crossprod(residuals(model))
#   if (missing(idata)) idata <- NULL
#   if (missing(idesign)) idesign <- NULL
#   if (!is.null(idata)){
#       for (i in 1:length(idata)){
#           if (is.null(attr(idata[,i], "contrasts"))){
#               contrasts(idata[,i]) <- if (is.ordered(idata[,i])) icontrasts[2]
#                   else icontrasts[1]
#           }
#       }
#       if (is.null(idesign)) stop("idesign (intra-subject design) missing.")
#       X.design <- model.matrix(idesign, data=idata)
#       if (check.imatrix) check.imatrix(X.design)
#       intercept <- has.intercept(X.design)
#       term.names <- term.names(idesign)
#       if (intercept) term.names <- c("(Intercept)", term.names)
#       which.terms <- match(iterms, term.names)
#       if (any(nas <- is.na(which.terms))){
#           if (sum(nas) == 1)
#               stop('The term "', iterms[nas],'" is not in the intrasubject design.')
#           else stop("The following terms are not in the intrasubject design: ",
#                   paste(iterms[nas], collapse=", "), ".")
#       }
#       select <- apply(outer(which.terms, attr(X.design, "assign") + intercept, "=="),
#           2, any)
#       P <- X.design[, select, drop=FALSE]
#   }
#   if (!is.null(P)){
#       rownames(P) <- colnames(B)
#       SSPE <- t(P) %*% SSPE %*% P
#       B <- B %*% P
#   }
#   rank <- sum(eigen(SSPE, only.values=TRUE)$values >= sqrt(.Machine$double.eps))
#   if (rank < ncol(SSPE))
#       stop("The error SSP matrix is apparently of deficient rank = ",
#           rank, " < ", ncol(SSPE))
#   r <- ncol(B)
#   if (is.null(rhs)) rhs <- matrix(0, nrow(L), r)
#   rownames(rhs) <- rownames(L)
#   colnames(rhs) <- colnames(B)
#   q <- NROW(L)
#   if (verbose){
#       cat("\nHypothesis matrix:\n")
#       print(L)
#       cat("\nRight-hand-side matrix:\n")
#       print(rhs)
#       cat("\nEstimated linear function (hypothesis.matrix %*% coef - rhs):\n")
#       print(drop(L %*% B - rhs))
#       cat("\n")
#   }
#   SSPH <- t(L %*% B - rhs) %*% solve(L %*% V %*% t(L)) %*% (L %*% B - rhs)
#   rval <- list(SSPH=SSPH, SSPE=SSPE, df=q, r=r, df.residual=df.residual, P=P,
#       title=title, test=test)
#   class(rval) <- "linearHypothesis.mlm"
#   rval
#}

print.linearHypothesis.mlm <- function(x, SSP=TRUE, SSPE=SSP,
		digits=getOption("digits"), ...){
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

linearHypothesis.survreg <- function(model, hypothesis.matrix, rhs=NULL,
		test=c("Chisq", "F"), vcov., verbose=FALSE, ...){
	if (missing(vcov.)) {
		vcov. <- vcov(model)
		p <- nrow(vcov.)
		vcov. <- vcov.[-p, -p]
	}
	linearHypothesis.default(model, hypothesis.matrix, rhs, test, vcov., verbose=verbose, ...)
}

linearHypothesis.polr <- function (model, hypothesis.matrix, rhs=NULL, vcov., verbose=FALSE, ...){
	k <- length(coef(model))
	V <- vcov(model)[1:k, 1:k]
	linearHypothesis.default(model, hypothesis.matrix, rhs, vcov.=V, verbose=verbose, ...)
}

coef.multinom <- function(object, ...){
	b <- nnet:::coef.multinom(object, ...)
	cn <- colnames(b)
	rn <- rownames(b)
	b <- as.vector(t(b))
	names(b) <- as.vector(outer(cn, rn, function(c, r) paste(r, c, sep=":")))
	b
}
