# generalized Durbin-Watson statistic (J. Fox)

# last modified 29 July 2001 by J. Fox

durbin.watson <- function(model, ...){
  UseMethod("durbin.watson")
  }

durbin.watson.lm <- function(model, max.lag=1, simulate=T, reps=1000, 
    method=c("resample","normal")){
    method<-match.arg(method)
    residuals<-residuals(model)
    if (any(is.na(residuals))) stop ('residuals include missing values')
    n<-length(residuals)
    r<-dw<-rep(0, max.lag)
    den<-sum(residuals^2)
    for (lag in 1:max.lag){
        dw[lag]<-(sum((residuals[(lag+1):n] - residuals[1:(n-lag)])^2))/den
        r[lag]<-(sum(residuals[(lag+1):n]*residuals[1:(n-lag)]))/den
        }
    if (!simulate){
        result<-list(r=r, dw=dw)
        class(result)<-"durbin.watson"
        result
        }
        else {
            S<-summary(model)$sigma
            X<-model.matrix(model)
            mu<-fitted.values(model)
            Y<-if (method == "resample") 
                matrix(sample(residuals, n*reps, replace=T), n, reps) + matrix(mu, n, reps)
                else matrix(rnorm(n*reps, 0, S), n, reps) + matrix(mu, n, reps)
            E<-residuals(lm(Y~X-1))
            DW<-apply(E, 2, durbin.watson, max.lag=max.lag)
            if (max.lag == 1) DW <- rbind(DW)
            p<-rep(0, max.lag)
            for (lag in 1:max.lag) {
                p[lag] <- (sum(dw[lag] < DW[lag,]))/reps
                p[lag] <- 2*(min(p[lag], 1-p[lag]))
                }
            result<-list(r=r, dw=dw, p=p)
            class(result)<-"durbin.watson"
            result
            }
    }

durbin.watson.default<-function(residuals, max.lag=1){
    if ( (!is.vector(residuals)) || (!is.numeric(residuals)) ) stop("requires vector of residuals")
    if (any(is.na(residuals))) stop ('residuals include missing values')
    n<-length(residuals)
    dw<-rep(0, max.lag)
    den<-sum(residuals^2)
    for (lag in 1:max.lag){
        dw[lag]<-(sum((residuals[(lag+1):n] - residuals[1:(n-lag)])^2))/den
        }
    dw
    }
    
print.durbin.watson<-function(x){
    max.lag<-length(x$dw)
    result<- if (is.null(x$p)) cbind(lag=1:max.lag,Autocorrelation=x$r, "D-W Statistic"=x$dw)
            else cbind(lag=1:max.lag,Autocorrelation=x$r, "D-W Statistic"=x$dw, "p-value"=x$p)
    rownames(result)<-rep("", max.lag)
    print(result)
    invisible(x)
    }
