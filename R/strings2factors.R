strings2factors <- function(object, which, not, exclude.unique, levels, verbose, ...){
  UseMethod("strings2factors")
}

strings2factors.data.frame <- function(object, which, not, exclude.unique=TRUE, levels=list(), 
                                        verbose=TRUE, ...){
  if (missing(which)) which <- sapply(object, is.character)
  if (is.numeric(which)) which <- names(object)[which]
  if (!missing(not) && !is.character(not)) not <- names(object)[not]
  if (!is.character(which)) which <- names(object)[which]
  if (!missing(not)) which <- setdiff(which, not)
  n <- nrow(object)
  for (var in which){
    levs <- levels[[var]]
    all.unique <- !anyDuplicated(object[[var]], incomparables=NA)
    if (all.unique){
      if (exclude.unique) {
        which <- setdiff(which, var)
        next
      }
      else warning("all values of ", var, " are unique")
    }
    object[[var]] <- if (is.null(levs)){
      factor(object[[var]])
    } else {
      factor(object[[var]], levels=levels[[var]])
    }
  }
  if (verbose){
    if (length(which) > 1){
      cat("\nThe following character variables were converted to factors\n", 
          paste(which, sep=", "), "\n")
    } else {
      cat("\n", which, "was converted to a factor")
    }
  }
  object
}
  
