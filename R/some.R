# adapted from head() and tail()
# 3/10/2017:  S. Weisberg modified to add an argument 'cols'
#             cols = num will display only the first num cols


some <- function(x, ...) UseMethod("some")

some.default <- function(x, n=10, ...){
    len <- length(x)
    ans <- x[sort(sample(len, min(n, len)))]
    if (length(dim(x)) == 1)
        array(ans, n, list(names(ans)))
    else ans
    }

some.matrix <- function(x, n=10, cols=NULL, ...){
  nr <- nrow(x)
  nc <- ncol(x)
  cols <- if(is.null(cols)) 1:nc else cols
  x[sort(sample(nr, min(n, nr))), cols, drop = FALSE]
  }

some.data.frame <- function(x, n=10, cols=NULL, ...){
    nr <- nrow(x)
    nc <- ncol(x)
    cols <- if(is.null(cols)) 1:nc else cols
    x[sort(sample(nr, min(n, nr))), cols, drop=FALSE]
    }
