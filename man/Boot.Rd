\name{Boot}
\alias{Boot}
\alias{Boot.default}
\alias{Boot.lm}
\alias{Boot.glm}
\alias{Boot.nls}

\title{Bootstrapping for regression models }
\description{
This function provides a simple front-end to the \code{boot} function in the
package also called \code{boot}.  Whereas \code{boot} is very general and therefore
has many arguments, the \code{Boot} function has very few arguments, but should
meet the needs of many users.
}
\usage{
\S3method{Boot}{default}(object, f=coef, labels=names(coef(object)),
                     R=999, method=c("case", "residual"))
                     
\S3method{Boot}{lm}(object, f=coef, labels=names(coef(object)),
                     R=999, method=c("case", "residual"))

\S3method{Boot}{glm}(object, f=coef, labels=names(coef(object)),
                     R=999, method=c("case", "residual"))

\S3method{Boot}{nls}(object, f=coef, labels=names(coef(object)),
                     R=999, method=c("case", "residual"))

}
\arguments{
  \item{object}{ A regression object of class \code{lm}, \code{glm} or \code{nls}.  The function may work with other regression objects that support the \code{update}
method and have a \code{subset} argument}
  \item{f}{A function whose one argument is the name of a regression object that will be applied to the updated regression object to
compute the statistics of interest.  The default is \code{coef}, to return to
regression coefficient estimates.  For example, 
\code{f = function(obj) coef(obj)[1]/coef(obj[2]} will bootstrap the ratio of the first and second coefficient estimates.}
  \item{labels}{Provides labels for the statistics computed by \code{f}.  If
this argument is of the wrong length, then generic labels will be generated.
}
  \item{R}{Number of bootstrap samples.  The number of bootstrap samples actually computed may be smaller than
this value if either the fitting method is iterative, or if the rank of a fittle \code{lm} or \code{glm}
model is different in the bootstrap replication than in the original data.}
  \item{method}{The bootstrap method, either \dQuote{case} for resampling
cases or \dQuote{residuals} for a residual bootstrap.  See the details
below. The residual bootstrap is available only for \code{lm} and \code{nls} objects and will return an error for \code{glm} objects.}
}
\details{
Whereas the \code{boot} function is
very general, \code{Boot} is very specific.  It takes the information from a
regression object and the choice of \code{method}, and creates a function that is
passed as the \code{statistic} argument to \code{boot}.  The argument \code{R}
is also passed to \code{boot}.  All other arguments to \code{boot} are
kept at their default values.

The methods available for \code{lm} and \code{nls} objects are \dQuote{case} and
\dQuote{residual}.  The case bootstrap resamples from the joint distribution
of the terms in the model and the response.  The residual bootstrap fixes the
fitted values from the original data, and creates bootstraps by adding a
bootstrap sample of the residuals to the fitted values to get a bootstrap
response.  It is an implementation of Algorithm 6.3, page 271, of
Davison and Hinkley (1997).  For \code{nls} objects ordinary residuals are used
in the resampling rather than the standardized residuals used in the \code{lm}
method.  The residual bootstrap for
generalized linear models has several competing approaches, but none are
without problems.  If you want to do a residual bootstrap for a glm, you
will need to write your own call to \code{boot}.

An attempt to model fit to a bootstrap sample may fail.  In a \code{lm} or
\code{glm} fit, the bootstrap sample could have a different rank from the original
fit.  In an \code{nls} fit, convergence may not be obtained for some bootstraps.
In either case, \code{NA} are returned for the value of the function \code{f}.
The summary methods handle the \code{NA}s appropriately.
}
\value{
See \code{\link[boot]{boot}} for the returned value from this function.  The car
package includes additional generic functions summary, confint and hist that works
with boot objects.
}
\references{
Davison, A, and Hinkley, D.  (1997)  \emph{Bootstrap Methods and their
Applications}.  Oxford:  Oxford University Press.

Fox, J. and Weisberg, S. (2011) \emph{Companion to Applied Regression}, Second Edition.
Thousand Oaks:  Sage.

Fox, J. and Weisberg, S. (2012) \emph{Bootstrapping}, 
\url{http://socserv.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Bootstrapping.pdf}.

Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley   
Wiley, Chapters 4 and 11.}

\author{ Sanford Weisberg, \email{sandy@umn.edu}.}

\seealso{Functions that work with
\code{Boot} objects from the \code{boot} package are 
\code{\link[boot]{boot.array}},
\code{\link[boot]{boot.ci}}, 
\code{\link[boot]{plot.boot}} and
\code{\link[boot]{empinf}}.  Additional
functions in the \code{car} package are 
\code{\link{summary.boot}},
\code{\link{confint.boot}}, and 
\code{\link{hist.boot}}.
}


\examples{
m1 <- lm(Fertility ~ ., swiss)
betahat.boot <- Boot(m1, R=199) # 199 bootstrap samples--too small to be useful
summary(betahat.boot)  # default summary
confint(betahat.boot)
hist(betahat.boot)
# Bootstrap for the estimated residual standard deviation:
sigmahat.boot <- Boot(m1, R=199, f=sigmaHat, labels="sigmaHat")
summary(sigmahat.boot)
confint(sigmahat.boot)
}
\keyword{regression}% at least one, from doc/KEYWORDS

