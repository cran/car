# Generalized Variance-Inflation Factors (J. Fox)

vif<-function(mod){
    #last modified 13 Dec 2000 by J. Fox
    UseMethod("vif")
    }

vif.lm<-function(mod) {
    #last modified 2 Dec 2003 by J. Fox
    if (!is.null(weights(mod))) stop("requires unweighted lm")
    if(!has.intercept(mod)) stop("requires model with intercept.")   
    terms<-term.names(mod)[-1]
    n.terms<-length(terms)
    if (n.terms < 2) stop("model contains fewer than 2 terms") 
    R<-cor(model.matrix(mod)[,-1])
    detR<-det(as.matrix(R))
    result<-matrix(0,n.terms,3)
    rownames(result)<-terms
    colnames(result)<-c("GVIF","Df","GVIF^(1/2Df)")
    assign<-mod$assign
    for (term in 1:n.terms){
        subs<-which(assign==term)-1
        result[term,1]<-det(as.matrix(R[subs,subs]))*
            det(as.matrix(R[-subs,-subs]))/detR
        result[term,2]<-length(subs)
        }
    if (all(result[,2]==1)) result<-result[,1]
        else result[,3]<-result[,1]^(1/(2*result[,2]))
    result
    }

vif.default<-function(mod){
    #last modified 13 Dec 2000 by J. Fox
    stop("requires lm object")
    }
