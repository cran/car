\name{Predict}
\alias{Predict}
\alias{Predict.lm}

\title{
Model Predictions
}

\description{
\code{Predict} is a generic function with, at present, a single method for \code{"lm"} objects,
\code{Predict.lm}, which is a modification of the standard \code{\link{predict.lm}} method in
the \pkg{stats} package, but with an additional \code{vcov.} argument for a user-specified covariance matrix for intreval estimation.
}

\usage{
Predict(object, ...)

\method{Predict}{lm}(object, newdata, se.fit = FALSE, 
  scale = NULL, df = Inf, 
  interval = c("none", "confidence", "prediction"), 
  level = 0.95, type = c("response", "terms"), 
  terms = NULL, na.action = na.pass,
  pred.var = res.var/weights, weights = 1, vcov., ...)
}

\arguments{
  \item{object}{a model object for which predictions are desired.}
  \item{newdata, se.fit, scale, df, interval, level, type, terms, na.action, pred.var, weights}{see \code{\link{predict.lm}}.}
  \item{vcov.}{optional, either a function to compute the coefficient covariance matrix of 
    \code{object} (e.g., \code{\link{hccm}}) or a coefficient covariance matrix (as returned,
    e.g., by \code{\link{hccm}}).}
  \item{\dots}{arguments to pass down to \code{Predict} or \code{predict} methods.}
}

\details{
If there is no appropriate method for \code{Predict}, then a \code{\link{predict}} method
is invoked. If there is a \emph{specific} \code{predict} method for the primary class of \code{object} but
only an \emph{inherited} \code{Predict} method, then the \code{predict} method is invoked.
Thus an object of class \code{c("glm", "lm")} will invoke method \code{predict.glm} rather than
\code{Predict.lm}, but an object of class \code{c("aov", "lm")} will invoke \code{Predict.lm}
rather than \code{predict.lm}.
}

\value{
See \code{\link{predict}} and \code{\link{predict.lm}}.
}

\references{
Fox, J. and Weisberg, S. (2019)
\emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{
John Fox \email{jfox@mcmaster.ca}
}

\seealso{
\code{\link{predict}}, \code{\link{predict.lm}}
}

\examples{
mod <- lm(interlocks ~ log(assets), data=Ornstein)
newd <- data.frame(assets=exp(4:12))
(p1 <- predict(mod, newd, interval="prediction"))
p2 <- Predict(mod, newd, interval="prediction", vcov.=vcov)
all.equal(p1, p2) # the same

(predict(mod, newd, se=TRUE))
(p3 <- Predict(mod, newd, se=TRUE, vcov.=hccm)) # larger SEs
p4 <- Predict(mod, newd, se=TRUE, vcov.=hccm(mod, type="hc3"))
all.equal(p3, p4) # the same
}

\keyword{models}
