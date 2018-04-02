# function Import March 14, 2017
# This add two arguments to the `import` file in the rio package
# import is just a front end to a number of file reading files and packages in R
# 3/14/2017:  S. Weisberg, wrote the file, that adds
#             row.names=TRUE, default, will select the left-most column of character data in the data file as
#             row names subject to length(x) == length(unique(x))
#             charAsFactor=TRUE converts character to factor if length(x) > length(unique(x))
#             logicalAsFactor=charAsFactor converts logical to factor
#             These arguments are read only if format %in% c("txt", "csv", "xls", "xlsx", "ods").
# 4/2/2017:  S. Weisberg changed and simplified arguments.
# 5/22/2017: S. Weisberg, fixed bug reading files with one character column (added drop=FALSE)

Import <- function(file, format, ..., row.names=TRUE,
                   stringsAsFactors = default.stringsAsFactors()){
  d <- rio::import(file, format, ...)
  fmt <- if(!missing(format)) format else{
    pos <- regexpr("\\.([[:alnum:]]+)$", file)
    ifelse(pos > -1L, substring(file, pos + 1L), "")
  }
# check for rows with no data
  d <- d[!apply(d, 1, function(row) all(is.na(row))), ]
  if(fmt %in% c("txt", "csv", "xls", "xlsx", "ods")){
    classes <- unlist(lapply(as.list(d), class))
    char <- classes %in% c("character", "logical")
    if(!any(char)) return(d)
    allUnique <- rep(FALSE, dim(d)[2])
      allUnique[char] <- apply(d[, char, drop=FALSE], 2, function(x) length(x) == length(unique(x)))
      if(row.names == TRUE){
        if(any(allUnique)){
          row.namesCol <- which(allUnique)[1] # use first non-repeating character col as row.names
          row.names(d) <- d[[row.namesCol]]   # set the row.names
          d <- d[, -row.namesCol]             # delete row.names column from data.frame
          allUnique <- allUnique[-row.namesCol]
          char <- char[-row.namesCol]}
      }
      if(stringsAsFactors & any(!allUnique)){
        for(j in which(char & !allUnique)) d[, j] <- factor(d[, j])
      }}
  return(d)
}
