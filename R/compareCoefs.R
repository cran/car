# 21 May 2010: small changes to output when there is just one model. J. Fox 
# 15 Aug 2010: changed name of function to compareCoefs to avoid name clash. J. Fox
# 18 May 2011: check for 'mer' objects, and handle them correctly. S. Weisberg
#  8 Sep 2011: check for 'lme' objects, and handle them correctly. S. Weisberg
# 11 Jan 2012: fix to work with any 'S4' object with a coef() method. 
#   suggested by David Hugh-Jones  University of Warwick http://davidhughjones.googlepages.com 
# 3 May 2012:  fixed bug if models are less than full rank.
# 17 Seot 2912: suppressing printing calls when there are none. J. Fox

compareCoefs <- function(..., se=TRUE, print=TRUE, digits=3){
    fixefmer <- function(m) {
      if(inherits(m, "mer")) m@fixef else fixef(m)
    }
    models <- list(...)
    n.models <- length(models)
    if (n.models < 1) return(NULL)
    getnames <- function(model) {
      if(inherits(model, "mer") | inherits(model, "lme")) names(fixef(model)) else
         names(coef(model))
         }
    getcoef <- function(model) {
      if(inherits(model, "mer") | inherits(model, "lme")) fixef(model) else 
          coef(model)
       }
    getcall <- function(model) {
      deparse(if (isS4(model)) model@call else model$call, 
               width.cutoff =  getOption("width") - 9)
       }
    getvar <- function(model) {
       if(inherits(model, "mer")) as.matrix(vcov(model)) else vcov(model)
       }
    coef.names <- unique(unlist(lapply(models, getnames)))
    table <- matrix(NA, length(coef.names), n.models*(1 + se))
    rownames(table) <- coef.names
    colnames(table) <- if (se) if (n.models > 1) paste(rep(c("Est.", "SE"), n.models),
                           rep(1:n.models, each=2)) else c("Estimate", "Std. Error")
        else if (n.models > 1) paste(rep("Est.", n.models), 1:n.models) else "Estimate" 
    calls <- !any(sapply(models, getcall) == "NULL")
    if(print == TRUE && calls) cat("\nCall:")
    for (i in 1:n.models){
        model <- models[[i]]
        fout <- deparse(getcall(model), width.cutoff=getOption("width") - 9)
		mod <- if (n.models > 1) paste(i, ":", sep="") else ""
        if(print == TRUE && calls) cat(paste("\n", mod, fout[1], sep=""))
        if(length(fout) > 1) for (f in fout[-1]) cat("\n",f)
        if (se) {
          ests <- getcoef(model)
          new <- cbind(ests, rep(NA, length(ests)))
          new[!is.na(ests), 2] <- sqrt(diag(getvar(model)))
          table[getnames(model), 2*(i - 1) + c(1, 2)] <- new}
#                cbind(getcoef(model), sqrt(diag(getvar(model)))) }
        else table[getnames(model), i] <- getcoef(model)
    }
    if(print == TRUE){ 
      cat("\n")
      printCoefmat(table, na.print="", digits=digits, tst.ind=NULL)} else
    table
}


