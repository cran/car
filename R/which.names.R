# positions of names in a data frame (J. Fox)

# last modified 9 Mar 2001 by J. Fox

which.names<-function(names, object){
    row.names<-if (inherits(object, "data.frame")) row.names(object) else object
    check<-outer(row.names, names, '==')
    if (!all(matched <- apply(check, 2, any))) 
        warning(paste(paste(names[!matched], collapse=", "), "not matched"))
    which(apply(check, 1, any))
    }
