Tapply <- function(formula, fun, data, na.action=na.pass, ..., targs=list()){
  yx <- if (missing(data)) model.frame(formula, na.action=na.action)
    else model.frame(formula, data=data, na.action=na.action)
  if (ncol(yx) < 2) stop("fewer than two variables")
  targs[c("X", "INDEX", "FUN")] <- list(yx[, 1], yx[, -1], fun)
  targs <- c(targs, list(...))
  do.call(tapply, targs)
}