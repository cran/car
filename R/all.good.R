# non-missing observations

# last modified 5 April 2001 by J. Fox

all.good <- function(...){
    apply(!is.na(cbind(...)), 1, all)
    }
