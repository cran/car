# Bonferroni test for an outlier (J. Fox)

outlier.test<-function(model, ...){
    UseMethod("outlier.test")
    }

outlier.test.lm<-function(model, labels=names(rstud)){
    #last modified 27 Jan 2001 by J. Fox
    rstud<-abs(rstudent(model))
    labels<-if(is.null(labels)) seq(along=rstud) else labels
    if (length(rstud) != length(labels)) 
        stop("Number of labels does not correspond to number of residuals.")
    rstud.max<-max(rstud)
    which.max<-which(rstud==rstud.max)
    df<-df.residual(model)-1
    n<-length(rstud)
    p<-2*(1-pt(rstud.max,df))
    result<-c(rstud.max, df, p, n*p)
    names(result)<-c("max|rstudent|", "df", "unadjusted p", "Bonferroni p")
    result<-list(test=result, obs=labels[which.max])
    class(result)<-"outlier.test"
    result
    }
    
outlier.test.glm<-function(model, labels=names(rstud)){
    #last modified 27 Jan 2001 by J. Fox
    rstud<-abs(rstudent(model))
    labels<-if(is.null(labels)) seq(along=rstud) else labels
    if (length(rstud) != length(labels)) 
        stop("Number of labels does not correspond to number of residuals.")
    rstud.max<-max(rstud)
    which.max<-which(rstud==rstud.max)
    n<-length(rstud)
    p<-2*(1-pnorm(rstud.max))
    result<-c(rstud.max, p, n*p)
    names(result)<-c("max|rstudent|", "unadjusted p", "Bonferroni p")
    result<-list(test=result, obs=labels[which.max])
    class(result)<-"outlier.test"
    result
    }

    
print.outlier.test<-function(x){
    # last modified 27 Jan 2001 by J. Fox
    test<-matrix(x$test, nrow=1)
    colnames(test)<-names(x$test)
    rownames(test)<-""
    print.matrix(test)
    if(length(x$obs)>1) cat("\nObservations:",x$obs,"\n")
      else cat("\nObservation:",x$obs,"\n")
    }
