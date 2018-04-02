# Export March 2017
# Wrapper for the 'export' function in the 'rio' package that adds automatic support for row names by
# converting the row names to the left-most column of the saved data file.
# 3/15/2017:  S. Weisberg

Export <- function(x, file, format, ..., keep.row.names){
  stopifnot(is.data.frame(x))
  if(!missing(keep.row.names)){
    name <- if(is.logical(keep.row.names) & keep.row.names==TRUE) "id" else keep.row.names
    if(name != FALSE){
      if(name %in% names(x))
        stop("There is a column named ", name, " already!")
      x <- cbind(rownames(x), x)
      names(x)[1] <- name
      attr(x, "row.names") <- 1:dim(x)[1]
  }}
  rio::export(x, file, format, ...)
}

