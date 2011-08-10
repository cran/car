# 21 May 2010: small changes to output when there is just one model. J. Fox 
# 15 Aug 2010: changed name of function to compareCoefs to avoid name clash. J. Fox
# 18 May 2011: check for 'mer' objects, and handle them correctly.

compareCoefs <- function(..., se=TRUE, digits=3){
    models <- list(...)
    n.models <- length(models)
    if (n.models < 1) return(NULL)
    getnames <- function(model) {
      if(inherits(model, "mer")) names(fixef(model)) else
         names(coef(model))
         }
    getcoef <- function(model) {
       if(inherits(model, "mer")) fixef(model) else coef(model)
       }
    getcall <- function(model) {
       if(inherits(model, "mer")) model@call else model$call
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
    cat("\nCall:")
    for (i in 1:n.models){
        model <- models[[i]]
        fout <- deparse(getcall(model), width.cutoff=getOption("width") - 9)
		mod <- if (n.models > 1) paste(i, ":", sep="") else ""
        cat(paste("\n", mod, fout[1], sep=""))
        if(length(fout) > 1) for (f in fout[-1]) cat("\n",f)
        if (se) {
          table[getnames(model), 2*(i - 1) + c(1, 2)] <-
                cbind(getcoef(model), sqrt(diag(getvar(model)))) }
        else table[getnames(model), i] <- getcoef(model)
    }
    cat("\n")
    printCoefmat(table, na.print="", digits=digits, tst.ind=NULL)
}
