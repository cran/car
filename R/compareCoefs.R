# 21 May 2010: small changes to output when there is just one model. J. Fox 
# 15 Aug 2010: changed name of function to compareCoefs to avoid name clash. J. Fox
# 18 May 2011: check for 'mer' objects, and handle them correctly. S. Weisberg
#  8 Sep 2011: check for 'lme' objects, and handle them correctly. S. Weisberg
# 11 Jan 2012: fix to work with any 'S4' object with a coef() method. 
#   suggested by David Hugh-Jones  University of Warwick http://davidhughjones.googlepages.com 
# 3 May 2012:  fixed bug if models are less than full rank.
# 17 Sept 2012: suppressing printing calls when there are none. J. Fox
# 22 June 2013: tweaks for lme4. J. Fox
# 26 Sept 2014: cleaned up printing of calls. J. Fox
# 2016-07-20: added test for model classes. J. Fox
# 2017-11-09: made calls to vcov() compatible with R 2.5.0. J. Fox
# 2017-12-05: added vcov. argument. J. Fox
# 2017-12-06,07: reorganized formatting of output, added zvals and pvals args. J. Fox
# 2018-02-07: removed leading blank line. J. Fox

compareCoefs <- function (..., se = TRUE, zvals = FALSE, pvals = FALSE, vcov., print = TRUE, digits = 3) {
  interleave <- function(X, se, zvals, pvals){
    krows <- 1 + se + zvals + pvals
    Y <- matrix(NA, nrow(X)*(krows + 1), ncol(X)/krows)
    nc <- ncol(Y)
    nr <- nrow(Y)
    rownames(Y) <- as.vector(rbind(rownames(X), 
                                   if (se) rep("SE", nrow(X)), 
                                   if (zvals) rep("z", nrow(X)),
                                   if (pvals) rep("Pr(>|z|)", nrow(X)),
                                   rep("", nrow(X))))
    colnames(Y) <- paste("Model", 1:nc)
    for (i in 1:nrow(X)){
      count <- 1
      Y[(i - 1)*(krows + 1) + count, ] <- X[i, (0:(nc - 1))*krows + count]
      if (se){
        count <- count + 1
        Y[(i - 1)*(krows + 1) + count, ] <- X[i, (0:(nc - 1))*krows + count]
      }
      if (zvals){
        count <- count + 1
        Y[(i - 1)*(krows + 1) + count, ] <- X[i, (0:(nc - 1))*krows + count]
      }
      if (pvals){
        count <- count + 1
        Y[(i - 1)*(krows + 1) + count, ] <- X[i, (0:(nc - 1))*krows + count]
      }
    }
    Y
  }
  splitExpr <- function(expr, width=getOption("width") - 4, at="[ ,=]"){
    if (length(grep("\n", expr)) >0 ){
      cmds <- strsplit(expr, "\n")[[1]]
      allcmds <- character(length(cmds))
      for (i in 1:length(cmds))
        allcmds[i] <- splitExpr(cmds[i], width=width, at=at)
      return(paste(allcmds, collapse="\n"))
    }
    if (nchar(expr) <= width) return(expr)
    where <- gregexpr(at, expr)[[1]]
    if (where[1] < 0) return(expr)
    singleQuotes <- gregexpr("'", expr)[[1]]
    doubleQuotes <- gregexpr('"', expr)[[1]]
    comment <- regexpr("#", expr)
    if (singleQuotes[1] > 0 && (singleQuotes[1] < doubleQuotes[1] || doubleQuotes[1] < 0 ) 
        && (singleQuotes[1] < comment[1] || comment[1] < 0 )){
      nquotes <- length(singleQuotes)
      if (nquotes < 2) stop("unbalanced quotes")
      for(i in seq(nquotes/2))
        where[(where > singleQuotes[2 * i - 1]) & (where < singleQuotes[2 * i])] <- NA
      where <- na.omit(where)
    }  
    else if (doubleQuotes[1] > 0 && (doubleQuotes[1] < singleQuotes[1] || singleQuotes[1] < 0) 
             && (doubleQuotes[1] < comment[1] || comment[1] < 0 )){
      nquotes <- length(doubleQuotes)
      if (nquotes < 2) stop("unbalanced quotes")
      for(i in seq(nquotes/2))
        where[(where > doubleQuotes[2 * i - 1]) & (where < doubleQuotes[2 * i])] <- NA
      where <- na.omit(where)
    }
    else if (comment > 0){
      where[where > comment] <- NA
      where <- na.omit(where)
    }
    if (length(where) == 0) return(expr)
    where2 <- where[where <= width]
    where2 <- if (length(where2) == 0) where[1]
    else where2[length(where2)]
    paste(substr(expr, 1, where2), "\n  ", 
          Recall(substr(expr, where2 + 1, nchar(expr)), width, at), sep="")
  } 
  removeExtraQuotes <- function(string) sub("\\\"$", "", sub("^\\\"", "", string))
  squeezeMultipleSpaces <- function(string) gsub(" {2,}", " ", string)
  intersection <- function(...){
    args <- list(...)
    if (length(args) == 2) intersect(args[[1]], args[[2]])
    else intersect(args[[1]], do.call(intersection, args[-1]))
  }
  
  models <- list(...)
  n.models <- length(models)
  if (n.models < 1) 
    return(NULL)
  if (n.models > 1){
    classes <- lapply(models, class)
    common.classes <- do.call(intersection, classes)
    if (length(common.classes) == 0) 
      warning("models to be compared are of different classes")
  }
  getnames <- function(model) {
    if (inherits(model, "merMod") || inherits(model, "mer") || 
        inherits(model, "lme")) 
      names(fixef(model))
    else names(coef(model))
  }
  getcoef <- function(model) {
    if (inherits(model, "merMod") || inherits(model, "mer") || 
        inherits(model, "lme")) 
      fixef(model)
    else coef(model)
  }
  getcall <- function(model) {
    paste(deparse(if (isS4(model))  model@call else model$call), collapse="")
  }
  getvar <- function(model) {
    if (inherits(model, "merMod") || inherits(model, "mer")) 
      as.matrix(vcov(model, complete=FALSE))
    else vcov(model, complete=FALSE)
  }
  vcov. <- if (missing(vcov.)) lapply(models, getvar)
  else{
    if (se) cat(paste(strwrap(paste0("\nStandard errors computed by ", 
                                     deparse(substitute(vcov.)), "\n\n"), 
                              width=getOption("width") - 2), collapse="\n  "))
    if (length(vcov.) == 1){
      if (!is.function(vcov.)) stop("vcov. is not a function")
      lapply(models, vcov.)
    }
    else{
      if (length(vcov.) != length(models))
        stop("number of entries in vcov. not equal to number of models")
      else {
        vc <- vector(length(models), mode="list")
        for (i in 1:length(models)){
          vc[[i]] <- if (is.function(vcov.[[i]])) vcov.[[i]](models[[i]]) 
          else if (is.matrix(vcov.[[i]])) vcov.[[i]] 
          else stop(i, "th element of vcov. is not a function or a matrix")
        }
        vc
      }
    }
  }
  coef.names <- unique(unlist(lapply(models, getnames)))
  table <- matrix(NA, length(coef.names), n.models * 4)
  rownames(table) <- coef.names
  colnames(table) <- if (n.models > 1) 
      paste(rep(c("Model", "SE", "z", "Pr(>|z|)"), n.models), rep(1:n.models, each = 4))
  else c("Estimate", "Std. Error", "z", "Pr(>|z|)")
  calls <- !any(sapply(models, getcall) == "NULL")
  if (print == TRUE && calls) 
    cat("Calls:")
  for (i in 1:n.models) {
    model <- models[[i]]
    fout <- getcall(model)
    mod <- if (n.models > 1) 
      paste(i, ": ", sep = "")
    else ""
    if (print && calls){
      cat(splitExpr(squeezeMultipleSpaces(paste("\n", mod, removeExtraQuotes(fout[1]), 
                                                sep = ""))))
    }
    if (print && calls && length(fout) > 1) {
      for (f in fout[-1]) cat("\n", splitExpr(squeezeMultipleSpaces(removeExtraQuotes(f))))
    }
    ests <- getcoef(model)
    aliased <- is.na(ests)
    new <- matrix(NA, length(ests), 4)
    new[, 1] <- ests
    new[!is.na(ests), 2] <- sqrt(diag(vcov.[[i]]))
    new[, 3] <- new[, 1]/new[, 2]
    new[, 4] <- 2*pnorm(abs(new[, 3]), lower.tail=FALSE)
    new[aliased, 1] <- -Inf
    table[getnames(model), 4 * (i - 1) + 1:4] <- new
  }
  table <- table[, c(TRUE, se, zvals, pvals)]
  if (print) {
    cat("\n\n")
    if (se || zvals || pvals) table <- interleave(table, se, zvals, pvals)
    cs.inds <- vector(n.models, mode="list")
    posn.coef <- 0
    posn.z <- rep(0, n.models)
    posn.p <- rep(0, n.models)
    for (i in 1:length(coef.names)){
      posn.coef <- posn.coef + 1
      posn.se <- posn.coef + se
      cs.inds[[i]] <- c(posn.coef:posn.se)
      if (zvals) posn.z[i] <- posn.se + 1
      if (pvals) posn.p[i] <- posn.se + zvals + 1
      posn.coef <- posn.se + zvals + pvals + any(c(se, zvals, pvals))
    }
    table.f <- formatCompareCoefs(table, cs.inds=cs.inds, digits = digits, 
                          tst.ind = if (zvals) posn.z else NULL, 
                          P.values= if (pvals) posn.p else NULL)
    print.default(table.f, quote = FALSE, right = TRUE)
  }
  invisible(table)
}

formatCompareCoefs <- function (x, digits, cs.inds, tst.ind, P.values = NULL){
  # this function adapted from stats::printCoefmat()
  x <- t(x)
  dig.tst <- max(1L, min(5L, digits - 1L))
  d <- dim(x)
  nc <- d[2L]
  xm <- data.matrix(x)
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  for (cs.ind in cs.inds){
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs != 0])) 
        floor(log10(range(acs[acs != 0], finite = TRUE)))
      else 0
      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits - digmin)), digits = digits)
    }
  }
  if (length(tst.ind)) 
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst), 
                            digits = digits)
  if (!is.null(P.values)){
    Cf[, P.values] <- format.pval(xm[, P.values], digits=dig.tst, eps=.Machine$double.eps)
  }
  Cf[Cf == "NA" | is.na(xm)] <- ""
  Cf[is.infinite(xm)] <- "aliased"
  t(Cf)
}
