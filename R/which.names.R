# positions of names in a data frame (J. Fox)

# last modified 2018-01-30 by J. Fox

whichNames <- function(names, object, ...){
  UseMethod("whichNames", object)
}

which.names <- function(names, object, ...){
  UseMethod("whichNames", object)
}

whichNames.data.frame <- function(names, object, ...){
    row.names <- row.names(object)
    check <- outer(row.names, names, '==')
    if (!all(matched <- apply(check, 2, any))) 
        warning(paste(paste(names[!matched], collapse=", "), "not matched"))
    result <- which(apply(check, 1, any))
    names(result) <- row.names[result]
    result[names[matched]]
}
	
whichNames.default <- function(names, object, ...){
  obj.names <- names(object)
  check <- outer(obj.names, names, '==')
  if (!all(matched <- apply(check, 2, any))) 
    warning(paste(paste(names[!matched], collapse=", "), "not matched"))
  result <- which(apply(check, 1, any))
  names(result) <- obj.names[result]
  result[names[matched]]
}