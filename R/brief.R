# added 2017-11-19 by J. Fox
# 2017-11-20: made S() methods brief() methods. J. Fox
# 2017-11-22: fixed brief.lm() and brief.glm() for models with aliased coeffs. J. Fox
# 2017-11-22: fixed bugs in brief.data.frame(), improved brief() and brief.list(). J. Fox
# 2017-12-15--21: tweaks to brief.data.frame. J. Fox
# 2017-12-19: added head, tail args to brief.data.frame()
# 2018-02-10: tweak brief.glm() output formatting
# 2018-12-26: Changed the argument for brief.lm from vcov.=vcov to just vcov. If arge is
#             missing set vcov. = vcov(orject, complete=FALSE) to match brief.glm

brief <- function(object, ...){
    g <- options("max.print"=.Machine$integer.max)
    on.exit(options(g))
    UseMethod("brief")
}

brief.matrix <- function(object, rows=if(nr <= 10) c(nr, 0) else c(3, 2), ...){
    nr <- nrow(object)
    brief.data.frame(object, rows, ...)
}

brief.data.frame <- function(object, rows=if(nr <= 10) c(nr, 0) else c(3, 2),
                          cols, head=FALSE, tail=FALSE, elided=TRUE, classes=inherits(object, "data.frame"), ...){
    pad <- function(x, right=TRUE){
        nch <- nchar(x)
        maxch <- max(nch)
        if (classes) maxch <- max(maxch, 3)
        if (right) paste0(x, strrep(" ", maxch - nch)) else paste0(strrep(" ", maxch - nch), x)
    }
    find.max.cols <- function(object, first, last, end=2){
        ncol <- ncol(object)
        nrow <- nrow(object)
        rows <- if (nrow > first + last) c(1:first, (nrow - last + 1):nrow) else 1:nrow
        nrows <- length(rows)
        object <- object[rows, , drop=FALSE]
        for(i in 1:(ncol - end)){
            res <- capture.output(
                if ((i + end) < ncol)
                    cbind(object[ , c(1:i, (ncol - end + 1):ncol), drop=FALSE], ". . .")
                else object[ , c(1:i, (ncol - end + 1):ncol), drop=FALSE])
            if (length(res) > nrows + 1) {
                i <- i - 1
                break
            }
        }
        if (i < 1){
            i <- 1
            end <- end - 1
        }
        c(i, end)
    }
    if (!isFALSE(head)){
        rows <- if (isTRUE(head)) c(6, 0) else c(head, 0)
    }
    if (!isFALSE(tail)){
        rows <- if (isTRUE(tail)) c(0, 6) else c(0, tail)
    }
    xx <- object
    dim <- dim(object)
    nr <- nrow(object)
    nc <- ncol(object)
    first <- rows[1]
    last <- rows[2]
    if (missing(cols)){
      cols <- find.max.cols(object, first, last)
    }
    first.c <- cols[1]
    last.c <- cols[2]
    if ((first.c + last.c) == 0 || (first + last) == 0) {
        stop("nothing to show")
        return(invisible(xx))
    }
    e.rows <- nr - (first + last)
    e.cols <- nc - (first.c + last.c)
    cat(dim[1], "x", dim[2], class(object)[1])
    if (elided && e.rows + e.cols > 0){
        cat(" (")
        if (e.rows > 0) cat(e.rows, "rows")
        if (e.rows > 0 && e.cols > 0) cat(" and ")
        if (e.cols > 0) cat(e.cols, "columns")
        cat (" omitted)")
    }
    cat("\n")
    if (length(elided) == 1) elided <- rep(elided, 2)
    force(classes)
    char <- is.character(object)
    rnms <- rownames(object)
    if (is.null(rnms)) {
        rnms <- paste0("[", 1:nr, ",]")
        rnames <- FALSE
    }
    else rnames <- TRUE
    nch <- nchar(rnms)
    mch <- max(nch[if (last == 0) 1:first else c(1:first, (nr - last + 1):nr)])
    rnms <- if (rnames) paste0(rnms, sapply(pmax(mch - nch, 0), function(x) paste0(rep(" ", x), collapse="")))
    else paste0(sapply(pmax(mch - nch, 0), function(x) paste0(rep(" ", x), collapse="")), rnms)
    rownames(object) <- rnms
    if (is.null(colnames(object))) {
        colnames(object) <- paste0("[,", 1:nc, "]")
    }
    object <- as.data.frame(object)
    if (nr - (first + last) > 0) object <- object[c(1:first, (nr - last + 1):nr), ]
    elided.cols <- FALSE
    if (nc - (first.c + last.c) > 0) {
        elided.cols <- TRUE
        object.left <- if (first.c > 0) cbind(object[, 1:first.c, drop=FALSE], rep("", nrow(object)))
                        else matrix(rep("", nrow(object)))
        object <- if (last.c > 0) cbind(object.left, object[, (nc - last.c + 1):nc, drop=FALSE])
                        else object.left
        colnames(object)[first.c + 1] <- ". . ."
    }
    col.classes <- paste0("[", substring(sapply(object, class), 1, 1), "]")
    for (j in 1:ncol(object)){
        if (is.numeric(object[, j])) {
            object[, j] <- format(object[, j])
            object[, j] <- pad(object[, j], right=FALSE)
        }
        else if (is.factor(object[, j])) {
            object[, j] <- droplevels(object[, j])
            levels(object[, j]) <- pad(levels(object[, j]))
        }
        else if (is.character(object[, j])) object[, j] <- pad(object[, j])
    }
    cnms <- colnames(object)
    object <- format(object)
    if (classes) {
        if (nc - (first.c + last.c) > 0) col.classes[first.c + 1] <- ""
        object <- rbind(col.classes , object)
        rownames(object)[1] <-""
        first <- first + 1
        nr <- nr + 1
    }
    if (first - classes > 0) {
        print(object[1:first, ], quote=char && !elided.cols)
        if (nr - (first + last) > 0) {
            cat(". . .")
            nch <- nchar(cnms)
            cnms <- sapply(nch, function(x) paste0(rep(" ", x), collapse=""))
            colnames(object) <- cnms
            if (last > 0) print(object[(first + 1):(first + last), ], quote=char && !elided.cols)
        }
    }
    else{
        object[1 + classes, ] <- rep("", ncol(object))
        rownames(object)[1 + classes] <- ". . ."
        print(object, quote=char && !elided.cols)
    }
    invisible(xx)
}

brief.function <- function(object, rows=c(5, 3), elided=TRUE, ...){
    first <- rows[1]
    last <- rows[2]
    fn <- format(object)
    nr <- length(fn)
    if (nr <= first + last) print(fn)
    else {
        cat(paste0(deparse(substitute(object)), " <- ", paste(fn[1:first], collapse="\n")))
        cat(paste0("\n\n. . . ", if (elided) paste0("(", nr - first - last, " lines omitted)"), "\n\n"))
        cat(paste(fn[(nr - last + 1):nr], collapse="\n"))
        cat("\n")
    }
    invisible(object)
}

brief.list <- function(object, rows=c(2, 1), elided=TRUE, ...){
    xx <- object
    first <- rows[1]
    last <- rows[2]
    nr <- length(object)
    if (nr <= first + last) print(object)
    else{
        cat(length(object),"element list")
        if (is.null(names(object))){
          names(object) <- 1:nr
        }
        for (i in 1:first){
            cat(paste0("\n[[", names(object[i]), "]]\n"))
            brief(object[[i]], elided=elided)
        }
        cat(paste0("\n. . . ", if (elided) paste0("(", nr - first - last, " list elements omitted)"), "\n"))
        for (i in (nr - last + 1):nr){
            cat(paste0("\n[[", names(object[i]), "]]\n"))
            brief(object[[i]], elided=elided)
        }
    }
    invisible(xx)
}

brief.vector <- function(object, rows=c(2, 1), elided=TRUE, ...){
    first <- rows[1]
    last <- rows[2]
    result <- capture.output(object)
    nr <- length(result)
    if (nr <= first + last) print(object)
    else{
       cat(length(object),"element", class(object)[1], "vector")
       cat("\n", paste0(result[1:first], "\n"))
       cat(paste0("\n. . . ", if (elided) paste0("(", nr - first - last, " lines omitted)"), "\n"))
       cat("\n", paste0(result[(nr - last + 1):nr]), "\n")
    }
    invisible(object)
}

# brief.vector() isn't a method and isn't exported

brief.integer <- brief.numeric <- brief.character <- brief.vector

brief.factor<- function(object, rows=c(2, 1), elided=TRUE, ...){
  first <- rows[1]
  last <- rows[2]
  result <- capture.output(object)
  levels <- result[length(result)]
  result <- result[-length(result)]
  nr <- length(result)
  if (nr <= first + last) print(object)
  else{
    cat(length(object),"element factor")
    cat("\n", paste0(result[1:first], "\n"))
    cat(paste0("\n. . . ", if (elided) paste0("(", nr - first - last, " lines omitted)"), "\n"))
    cat("\n", paste0(result[(nr - last + 1):nr]), "\n")
    cat(levels)
  }
  invisible(object)
}

# methods for statistical models

brief.default <- function(object, terms = ~ ., intercept=missing(terms), pvalues=FALSE, digits=3, horizontal=TRUE, ...){
  sumry <- summary(object)
  if (is.atomic(object) || is.atomic(sumry) || is.null(sumry$coefficients) || !is.matrix(sumry$coefficients)){
    if (is.vector(object)) brief.vector(object, ...)
    else if (is.list(object)) brief.list(object, ...)
    else stop("no appropriate method for object of class '", class(object), "'")
    return(invisible(object))
  }
  use <- coefs2use(object, terms, intercept)
  cols <- if (pvalues) c(1, 2, 4) else 1:2
  coefs <- sumry$coefficients[use, cols, drop=FALSE]
  colnames(coefs) <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|z|)") else c("Estimate", "Std. Error")
  print(if (horizontal) t(coefs) else coefs, digits=digits)
  invisible(sumry)
}

brief.lm <- function(object, terms = ~ ., intercept=missing(terms), pvalues=FALSE, digits=3, horizontal=TRUE, vcov., ...){ 
  use <- coefs2use(object, terms, intercept)
  vcov. <- if(missing(vcov.)) vcov(object, complete=FALSE) else vcov.
  sumry <- S(object, vcov.=vcov., ...)
  cols <- if (pvalues) c(1, 2, 4) else 1:2
  coefs <- sumry$coefficients
  if (!is.null(aliased <- sumry$aliased) && any(aliased)) {
    cn <- names(aliased)
    coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
                                                            colnames(coefs)))
    coefs[!aliased, ] <- sumry$coefficients
  }
  coefs <- coefs[use, cols, drop=FALSE]
  n.aliased <- sum(is.na(coefs[, 1]))
  if (n.aliased > 0)  cat(n.aliased, if(n.aliased > 1) "coefficients" else "coefficient",
                          "not defined because of singularities\n\n")
  colnames(coefs) <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|t|)") else c("Estimate", "Std. Error")
  print(if (horizontal) t(coefs) else coefs, digits=digits)
  if (missing(terms)) cat("\n Residual SD =", format(sumry$sigma, digits=digits),
                          "on", object$df.residual, "df, R-squared =", format(sumry$r.squared, digits=digits))
  invisible(sumry)
}

brief.glm <- function(object, terms = ~ ., intercept=missing(terms), pvalues=FALSE, digits=3, horizontal=TRUE, vcov., dispersion, exponentiate, ...){
  if (!missing(vcov.) && !missing(dispersion))
    stop("cannot specify both the dispersion and vcov. arguments")
  if (missing(exponentiate)) exponentiate <- object$family$link %in% c("log", "logit")
  use <- coefs2use(object, terms, intercept)
  sumry <- if (!missing(vcov.)) S(object, digits, vcov.=vcov., ...)
  else if (!missing(dispersion)) S(object, digits, dispersion=dispersion, ...)
  else summary(object, ...)
  cols <- if (pvalues) c(1, 2, 4) else 1:2
  coefs <- sumry$coefficients
  if (!is.null(aliased <- sumry$aliased) && any(aliased)) {
    cn <- names(aliased)
    coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn,
                                                            colnames(coefs)))
    coefs[!aliased, ] <- sumry$coefficients
  }
  coefs <- coefs[use, cols, drop=FALSE]
  n.aliased <- sum(is.na(coefs[, 1]))
  if (n.aliased > 0)  cat(n.aliased, if(n.aliased > 1) "coefficients" else "coefficient",
                          "not defined because of singularities\n\n")
  colnames(coefs) <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|t|)") else c("Estimate", "Std. Error")
  if (exponentiate){
    coefs <- cbind(coefs, exp(coefs[, 1]))
    colnames(coefs)[if (pvalues) 4 else 3] <- "exp(Estimate)"
  }
  print(if (horizontal) t(coefs) else coefs, digits=digits)
  if (missing(terms)) cat(paste0("\n Residual deviance = ", format(object$deviance, digits=digits),
                          " on ", object$df.residual, " df",
                          if (family(object)$family %in% c("binomial", "poisson")) ""
                          else (paste0(", Est. dispersion = ", format(sumry$dispersion, digits=digits)))))
  invisible(sumry)
}

brief.polr <- function(object, terms = ~ ., intercept, pvalues=FALSE, digits=3, horizontal=TRUE, exponentiate=TRUE, ...){
  sumry <- summary(object)
  coefs <- sumry$coefficients[ , 1:2]
  if (pvalues) {
    coefs <- cbind(coefs, 2*pnorm(abs(coefs[ , 1]/coefs[, 2]), lower.tail=FALSE))
  }
  use <- if (missing(terms)) 1:nrow(coefs) else coefs2use(object, terms, FALSE)
  coefs <- coefs[use, , drop=FALSE]
  colnames(coefs) <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|z|)") else c("Estimate", "Std. Error")
  if (exponentiate){
    coefs <- cbind(coefs, exp(coefs[, 1]))
    colnames(coefs)[if (pvalues) 4 else 3] <- "exp(Estimate)"
    if (missing(terms)){
      n.thresholds <- length(object$zeta)
      n.pars <- nrow(coefs)
      coefs[(n.pars - n.thresholds + 1):n.pars , if (pvalues) 4 else 3] <- NA
    }
  }
  print(if (horizontal) t(coefs) else coefs, digits=digits, na.print="")
  if (missing(terms)) cat("\n Residual deviance =", format(object$deviance, digits=digits),
                          "on", object$df.residual, "df")
  invisible(sumry)
}

brief.multinom <- function(object, terms = ~ ., intercept=missing(terms), pvalues=FALSE, digits=3, horizontal=TRUE, exponentiate=TRUE, ...){
  use <- coefs2use(object, terms, intercept)
  sumry <- summary(object, ...)
  b <- sumry$coefficients
  se <- sumry$standard.errors
  p <- 2*pnorm(abs(b/se), lower.tail=FALSE)
  levels <- sumry$lev
  labels <- if (pvalues) c("Estimate", "Std. Error", "Pr(>|z|)") else c("Estimate", "Std. Error")
  if (exponentiate) labels <- c(labels, "exp(Estimate)")
  if (length(levels) == 2){
    b <- b[use]
    se <- se[use]
    p <- p[use]
    table <- if (pvalues) rbind(b, se, p) else rbind(b, se)
    if (exponentiate) table <- rbind(table, exp(b))
    rownames(table) <- labels
    cat("\n ", levels[2], "\n")
    print(if (horizontal) table else t(table), digits=digits)
  }
  else{
    b <- b[, use, drop=FALSE]
    se <- se[, use, drop=FALSE]
    p <- p[, use, drop=FALSE]
    table <- if (pvalues) abind(t(b), t(se), t(p), along=1.5) else abind(t(b), t(se), along=1.5)
    if (exponentiate) table <- abind(table, t(exp(b)), along=2)
    dimnames(table)[[2]] <- labels
    for (level in levels[-1]){
      cat("\n ", level, "\n")
      result <- if (horizontal) t(table[, , level]) else table[, , level]
      if (dim(table)[1] == 1){
        if (horizontal) rownames(result) <- dimnames(table)[1] else {
          result <- matrix(result, ncol=1)
          colnames(result) <- dimnames(table)[1]
        }
      }
      print(result, digits=digits)
    }
  }
  if (missing(terms)) cat("\n Residual deviance =", format(object$deviance, digits=digits),
                          "fitting", length(b), "parameters")
  invisible(sumry)
}

