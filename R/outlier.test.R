# Bonferroni test for an outlier (J. Fox)

# last modified 29 Jan 04 by J. Fox

outlier.test<-function(model, ...){
    UseMethod("outlier.test")
    }

outlier.test.lm<-function(model, labels=names(rstud), ...){
    #last modified 13 Nov 2001 by J. Fox
    rstud<-abs(rstudent(model))
    labels<-if(is.null(labels)) seq(along=rstud) else labels
    if (length(rstud) != length(labels)) 
        stop("Number of labels does not correspond to number of residuals.")
    rstud.max<-max(rstud, na.rm=TRUE)
    which.max<-which(rstud==rstud.max)
    df<-df.residual(model)-1
    n<-sum(!is.na(rstud))
    p<-2*(1-pt(rstud.max,df))
    bp <- if (n*p <= 1) n*p else NA
    result<-c(rstud.max, df, p, bp)
    names(result)<-c("max|rstudent|", "df", "unadjusted p", "Bonferroni p")
    result<-list(test=result, obs=labels[which.max])
    class(result)<-"outlier.test"
    result
    }
    
outlier.test.glm<-function(model, labels=names(rstud), ...){
    #last modified 13 Nov 2001 by J. Fox
    rstud<-abs(rstudent(model))
    labels<-if(is.null(labels)) seq(along=rstud) else labels
    if (length(rstud) != length(labels)) 
        stop("Number of labels does not correspond to number of residuals.")
    rstud.max<-max(rstud, na.rm=TRUE)
    which.max<-which(rstud==rstud.max)
    n<-sum(!is.na(rstud))
    p<-2*(1-pnorm(rstud.max))
    bp <- if (n*p <= 1) n*p else NA
    result<-c(rstud.max, p, bp)
    names(result)<-c("max|rstudent|", "unadjusted p", "Bonferroni p")
    result<-list(test=result, obs=labels[which.max])
    class(result)<-"outlier.test"
    result
    }

    
print.outlier.test<-function(x, digits=options("digits")[[1]], ...){
    # last modified 29 Jan 2004 by J. Fox
    test<-signif(x$test, digits=digits)
    if (length(test) == 4){
        cat(paste("\nmax|rstudent| = ", test[1], ", degrees of freedom = ", test[2],
            ",\nunadjusted p = ", test[3], 
            ", Bonferroni p", if (is.na(test[4])) " > 1" else paste(" =", test[4]), "\n",
            sep=""))
        }
    else {
        cat(paste("\nmax|rstudent| = ", test[1],
            ",\nunadjusted p = ", test[2], 
            ", Bonferroni p", if (is.na(test[3])) " > 1" else paste(" =", test[3]), "\n",
            sep=""))
        }        
    if(length(x$obs)>1) cat("\nObservations:",x$obs,"\n")
      else cat("\nObservation:",x$obs,"\n")
    invisible(x)
    }
