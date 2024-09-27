# added 2022-05-24 by J. Fox
# 2024-04-11: invisibly return coordinates. J. Fox

avPlot3d <- function(model, coef1, coef2, id=TRUE, ...) {
  UseMethod("avPlot3d")
}

avPlot3d.lm <- function(model, coef1, coef2, id=TRUE, fit="linear", ...){
  fit <- match.arg(fit, c("linear", "robust"), several.ok=TRUE)
  coefs <- names(coef(model))
  which1 <- which(coef1 == coefs)
  if (length(which1) == 0) stop(coef1, " is not in the model")
  which2 <- which(coef2 == coefs)
  if (length(which2) == 0) stop(coef2, " is not in the model")
  y <- responseName(model) 
  X <- model.matrix(model)
  D <- data.frame(
    x1 = residuals(lm(X[, which1] ~ X[, -c(which1, which2)] - 1)),
    x2 = residuals(lm(X[, which2] ~ X[, -c(which1, which2)] - 1)),
    y1 = residuals(lm(model.response(model.frame(model)) ~ X[, -c(which1, which2)] - 1))
  )
  scatter3d(y1 ~ x1 + x2, data=D, xlab=paste0(coef1, " | others"), 
            zlab=paste0(coef2, " | others"), ylab=paste0(y, " | others"), id=id, fit=fit, ...)
  colnames(D) <- c(coef1, coef2, y)
  invisible(D)
}

avPlot3d.glm <- function(model, coef1, coef2, id=TRUE, type=c("Wang", "Weisberg"), 
                         fit="linear", ...){
  type <- match.arg(type)
  fit <- match.arg(fit, c("linear", "robust"), several.ok=TRUE)
  coefs <- names(coef(model))
  which1 <- which(coef1 == coefs)
  if (length(which1) == 0) stop(coef1, " is not in the model")
  which2 <- which(coef2 == coefs)
  if (length(which2) == 0) stop(coef2, " is not in the model")
  y <- responseName(model) 
  X <- model.matrix(model)
  wt <- model$prior.weights
  wt2 <- if (type == "Wang") wt*model$weights else wt
  D <- data.frame(
    x1 = residuals(lm(X[, which1] ~ X[, -c(which1, which2)] - 1, weights=wt2)),
    x2 = residuals(lm(X[, which2] ~ X[, -c(which1, which2)] - 1, weights=wt2)),
    y1 = residuals(glm(model.response(model.frame(model)) ~ X[, -c(which1, which2)] - 1,
                       family=family(model), weights=wt), 
                   type="pearson")
  )
  scatter3d(y1 ~ x1 + x2, data=D, xlab=paste0(coef1, " | others"), 
            zlab=paste0(coef2, " | others"), ylab=paste0(y, " | others"), id=id, fit=fit, ...)
  colnames(D) <- c(coef1, coef2, y)
  invisible(D)
}
