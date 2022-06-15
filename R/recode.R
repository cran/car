# recode function (J. Fox)
# 2019-11-14: change class(x) == "y" to inherits(x, "y")
# 2022-05-24: add to.value, interval, separator args

recode <- function(var, recodes, as.factor, as.numeric=TRUE, levels, 
                   to.value="=", interval=":", separator=";"){
  squeezeBlanks <- function(text){
    gsub(" *", "",  text)
  }
  if (grepl(paste0("[", separator, interval, "]"), to.value)) 
    stop("to.value may not contain ", interval, " or ", separator)
  lo <- -Inf
  hi <- Inf
  recodes <- gsub("\n|\t", " ", recodes)
  recode.list <- rev(strsplit(recodes, separator)[[1]])
  is.fac <- is.factor(var)
  if (missing(as.factor)) as.factor <- is.fac
  if (is.fac) var <- as.character(var)
  result <- var
  for (term in recode.list){
    if (0 < length(grep(interval, term))) {
      range <- strsplit(strsplit(term, to.value)[[1]][1],interval)
      low <- try(eval(parse(text=range[[1]][1])), silent=TRUE)
      if (inherits(low, "try-error")){
        stop("\n  in recode term: ", term, 
             "\n  message: ", low)
      }
      high <- try(eval(parse(text=range[[1]][2])), silent=TRUE)
      if (inherits(high, "try-error")){
        stop("\n  in recode term: ", term, 
             "\n  message: ", high)
      }
      target <- try(eval(parse(text=strsplit(term, to.value)[[1]][2])), silent=TRUE)
      if (inherits(target, "try-error")){
        stop("\n  in recode term: ", term, 
             "\n  message: ", target)
      }
      result[(var >= low) & (var <= high)] <- target
    }
    else if (0 < length(grep(paste0("^else", to.value), squeezeBlanks(term)))) {
      target <- try(eval(parse(text=strsplit(term, to.value)[[1]][2])), silent=TRUE)
      if (inherits(target, "try-error")){
        stop("\n  in recode term: ", term, 
             "\n  message: ", target)
      }
      result[1:length(var)] <- target
    }
    else {
      set <- try(eval(parse(text=strsplit(term, to.value)[[1]][1])), silent=TRUE)
      if (inherits(set, "try-error")){
        stop("\n  in recode term: ", term, 
             "\n  message: ", set)
      }
      target <- try(eval(parse(text=strsplit(term, to.value)[[1]][2])), silent=TRUE)
      if (inherits(target, "try-error")){
        stop("\n  in recode term: ", term, 
             "\n  message: ", target)
      }
      for (val in set){
        if (is.na(val)) result[is.na(var)] <- target
        else result[var == val] <- target
      }
    }
  }
  if (as.factor) {
    result <- if (!missing(levels)) factor(result, levels=levels) 
    else as.factor(result)
  }
  else if (as.numeric && (!is.numeric(result))) {
    result.valid <- na.omit(result)
    opt <- options("warn" = -1)
    result.valid <- as.numeric(result.valid)
    options(opt)
    if (!any(is.na(result.valid))) result <- as.numeric(result)
  }
  result
}

Recode <- function (...) car::recode(...)

