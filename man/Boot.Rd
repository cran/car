\name{Boot}
\alias{Boot}
\alias{Boot.default}
\alias{Boot.lm}
\alias{Boot.glm}
\alias{Boot.nls}

\title{Bootstrapping for regression models }
\description{
This function provides a simple front-end to the \code{boot} function in the
\pkg{boot} package that is tailored to bootstrapping based on regression models.  Whereas \code{boot} is very general and therefore
has many arguments, the \code{Boot} function has very few arguments.
}
\usage{
Boot(object, f=coef, labels=names(f(object)), R=999,
  method=c("case", "residual"), ncores=1, \dots)

\S3method{Boot}{default}(object, f=coef, labels=names(f(object)),
  R=999, method=c("case", "residual"), ncores=1,
  start = FALSE, \dots)

\S3method{Boot}{lm}(object, f=coef, labels=names(f(object)),
  R=999, method=c("case", "residual"), ncores=1, \dots)

\S3method{Boot}{glm}(object, f=coef, labels=names(f(object)),
  R=999, method=c("case", "residual"), ncores=1, \dots)

\S3method{Boot}{nls}(object, f=coef, labels=names(f(object)),
  R=999, method=c("case", "residual"), ncores=1, \dots)

}
\arguments{
\item{object}{A regression object of class \code{"lm"}, \code{"glm"} or 
   \code{"nls"}.  The function may work with other regression objects that support
   the \code{update} method and have a \code{subset} argument.  See discussion of
   the default method in the details below.}
\item{f}{A function whose one argument is the name of a regression object that will
   be applied to the updated regression object to compute the statistics of 
   interest.  The default is \code{coef}, to return regression coefficient 
   estimates.  For example, \code{f = function(obj) coef(obj)[1]/coef(obj)[2]} will
   bootstrap the ratio of the first and second coefficient estimates.}
\item{labels}{Provides labels for the statistics computed by \code{f}.  Default
   labels are obtained from a call to \code{f}, or generic labels if \code{f} 
   does not return names.}
\item{R}{Number of bootstrap samples.  The number of bootstrap samples actually
   computed may be smaller than this value if either the fitting method is 
   iterative and fails to converge for some boothstrap samples, or if the rank of a
   fitted model is different in a bootstrap replication than in the original data.}
\item{method}{The bootstrap method, either \dQuote{case} for resampling
   cases or \dQuote{residuals} for a residual bootstrap.  See the details
   below. The residual bootstrap is available only for \code{lm} and \code{nls}
   objects and will return an error for \code{glm} objects.}
\item{\dots}{Arguments passed to the \code{boot} function, see 
   \code{\link[boot]{boot}}.}
\item{start}{Should the estimates returned by \code{f} be passed as starting
    values for each bootstrap iteration? Alternatively, \code{start} can be a
    numeric vector of starting values.  The default is to use
    the estimates from the last bootstrap iteration as starting values
    for the next iteration.}
\item{ncores}{A numeric argument that specifies the number of cores for parallel
    processing for unix systems.  If less than or equal to 1, no parallel processing wiill be used. Note in a Windows platform will produce a warning and set this argument to 1.}
}
\details{
\code{Boot} uses a
regression object and the choice of \code{method}, and creates a function that is
passed as the \code{statistic} argument to the \code{boot} function in the \pkg{boot} package.  The argument \code{R} is also passed to \code{boot}.  If \code{ncores} is greater than 1, then the \code{parallel} and \code{ncpus} arguments to \code{boot} are set appropriately to use multiple codes, if available, on your computer. All other arguments to \code{boot} are kept at their default values unless you pass values for them.

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

For the default object to work with other types of regression models, the model must have methods for the the following generic functions:  \code{residuals(object, type="pearson")} must return Pearson residuals; \code{fitted(object)} must return fitted values; \code{hatvalues(object)} should return the leverages, or perhaps the value 1 which will effectively ignore setting the hatvalues.  In addition, the \code{data} argument should contain no missing values among the columns actually used in fitting the model, as the resampling may incorrectly attempt to include cases with missing values.  For \code{lm}, \code{glm} and \code{nls}, missing values cause the return of an error message.

An attempt to fit using a bootstrap sample may fail.  In a \code{lm} or
\code{glm} fit, the bootstrap sample could have a different rank from the original
fit.  In an \code{nls} fit, convergence may not be obtained for some bootstraps.
In either case, \code{NA} are returned for the value of the function \code{f}.
The summary methods handle the \code{NA}s appropriately.

Fox and Weisberg (2017) cited below discusses this function and provides more examples.
}
\value{
See \code{\link[boot]{boot}} for the returned value of the structure returned by this function.  
}
\references{
Davison, A, and Hinkley, D.  (1997)  \emph{Bootstrap Methods and their
Applications}.  Oxford:  Oxford University Press.

Fox, J. and Weisberg, S. (2019) \emph{Companion to Applied Regression}, Third Edition. Thousand Oaks:  Sage.

Fox, J. and Weisberg, S. (2019) \emph{Bootstrapping Regression Models in R},
\url{https://socialsciences.mcmaster.ca/jfox/Books/Companion/appendices/Appendix-Bootstrapping.pdf}.

Weisberg, S. (2014) \emph{Applied Linear Regression}, Fourth Edition, Wiley
Wiley, Chapters 4 and 11.}

\author{ Sanford Weisberg, \email{sandy@umn.edu}.  Achim Zeileis added multicore support, and also fixed the default method to work for many more regression models.}

\seealso{Functions that work with
\code{boot} objects from the \pkg{boot} package are
\code{\link[boot]{boot.array}},
\code{\link[boot]{boot.ci}},
\code{\link[boot]{plot.boot}} and
\code{\link[boot]{empinf}}.  Additional
functions in the \pkg{car} package are
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

