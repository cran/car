symbox <- function(x, powers = c(-1, -0.5, 0, 0.5, 1), start=0) {
    if (!(is.vector(x) && is.numeric(x))) stop("x should be a numeric vector.")
    x <- x + start
    if (min(x) <= 0) stop("Negative values for x are not allowed (after start applied).")
    result <- outer(x, powers,  "^")
    result[,powers < 0] <- - result[,powers < 0]
    result[,powers == 0] <- log(x)
    result <- as.data.frame(scale(result))
    names <- as.character(powers)
    names[powers == 0] <- "log"
    names(result) <- names
    boxplot(result, xlab="Powers")
}
