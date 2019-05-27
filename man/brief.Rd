\name{brief}
\alias{brief}
\alias{brief.data.frame}
\alias{brief.matrix}
\alias{brief.numeric}
\alias{brief.integer}
\alias{brief.character}
\alias{brief.factor}
\alias{brief.function}
\alias{brief.list}
\alias{brief.lm}
\alias{brief.glm}
\alias{brief.multinom}
\alias{brief.polr}
\alias{brief.default}

\title{
Print Abbreviated Ouput
}
\description{
Print data objects and statistical model summaries in abbreviated form.
}
\usage{
brief(object, ...)

\method{brief}{data.frame}(object, rows = if (nr <= 10) c(nr, 0) else c(3, 2),
    cols, head=FALSE, tail=FALSE, elided = TRUE,
    classes = inherits(object, "data.frame"), ...)
\method{brief}{matrix}(object, rows = if (nr <= 10) c(nr, 0) else c(3, 2), ...)

\method{brief}{numeric}(object, rows = c(2, 1), elided = TRUE, ...)
\method{brief}{integer}(object, rows = c(2, 1), elided = TRUE, ...)
\method{brief}{character}(object, rows = c(2, 1), elided = TRUE, ...)
\method{brief}{factor}(object, rows=c(2, 1), elided=TRUE, ...)

\method{brief}{list}(object, rows = c(2, 1), elided = TRUE, ...)

\method{brief}{function}(object, rows = c(5, 3), elided = TRUE, ...)

\method{brief}{lm}(object,  terms = ~ .,
    intercept=missing(terms), pvalues=FALSE,
    digits=3, horizontal=TRUE, vcov., ...)

\method{brief}{glm}(object, terms = ~ .,
    intercept=missing(terms), pvalues=FALSE,
    digits=3, horizontal=TRUE, vcov., dispersion, exponentiate, ...)

\method{brief}{multinom}(object, terms = ~ .,
    intercept=missing(terms), pvalues=FALSE,
    digits=3, horizontal=TRUE, exponentiate=TRUE, ...)

\method{brief}{polr}(object, terms = ~ .,
    intercept, pvalues=FALSE,
    digits=3, horizontal=TRUE, exponentiate=TRUE, ...)

\method{brief}{default}(object, terms = ~ .,
    intercept=missing(terms), pvalues=FALSE,
    digits=3, horizontal=TRUE, ...)
}

\arguments{
  \item{object}{a data or model object to abbreviate.}

  \item{rows}{for a matrix or data frame, a 2-element integer vector with the number of rows to print at
  the beginning and end of the display; for a vector or factor,
  the number of lines of output to show at the beginning
  and end; for a list, the number of elements to show at the beginning and end; for a function, the number
  of lines to show at the beginning and end.}

  \item{cols}{for a matrix or data frame, a 2-element integer vector with the number of columns to print at
  the beginning (i.e., left) and end (right) of the display.}

  \item{head, tail}{alternatives to the \code{rows} argument; if \code{TRUE}, print the first or last 6
  rows; can also be the number of the first or last few rows to print; only one of \code{heads} and
  \code{tails} should be specified; ignored if \code{FALSE} (the default).}

  \item{elided}{controls whether to report the number of elided elements, rows, or columns; default is \code{TRUE}.}

  \item{classes}{show the class of each column of a data frame at the top of the column; the classes are
  shown in single-character abbreviated form---e.g., \code{[f]} for a factor, \code{[i]} for an integer
  variable, \code{[n]} for a numeric variable, \code{[c]} for a character variable.}

    \item{terms}{a one-sided formula giving the terms to summarize; the default is \code{~ .}---i.e., to summarize all terms in the model.}
  \item{intercept}{whether or not to include the intercept; the default is \code{TRUE} unless the \code{terms} argument is given, in which
    case the default is \code{FALSE}; ignored for \code{polr} models.}
  \item{pvalues}{include the p-value for each coefficient in the table; default is \code{FALSE}.}
    \item{exponentiate}{for a \code{"glm"} or \code{"glmerMod"} model using the \code{log} or \code{logit} link, or a
    \code{"polr"} or \code{"multinom"} model, show exponentiated coefficient estimates and confidence bounds.}
    \item{digits}{significant digits for printing.}
     \item{horizontal}{if \code{TRUE} (the default), orient the summary produced by \code{brief} horizontally, which typically saves space.}
       \item{dispersion}{see \code{\link{summary.glm}}}
  
  \item{vcov.}{either a matrix giving the estimated covariance matrix of the estimates, 
    or a function that
    when called with \code{object} as an argument returns an estimated covariance matrix 
    of the estimates.  The default \code{vcov. = vcov(object, complete=FALSE)} uses the 
    usual estimated covariance matrix with NA for any row and column with aliased regressors. 
    Other choices include the functions documented at \code{\link{hccm}}, and a bootstrap
    estimate \code{vcov.=vcov(Boot(object))}; see the documentation for \code{\link{Boot}}.  
    For the \code{glm} method, if the \code{vcov.} or \code{dispersion} argument is specified,
    then Wald-based confidence limits are computed; otherwise the reported confidence limits
    are computed by profiling the likelihood.
    NOTES:  (1) The \code{dispersion} and \code{vcov.} arguments may not \emph{both} be 
    specified. (2) Setting \code{vcov.=vcov} returns an error if the model includes aliased 
    terms; use \code{vcov.=vcov(object, complete=FALSE)}. (3) The \code{hccm} method will 
    generally return a matrix of full rank even if the model has aliased terms.  Similarly
    \code{vcov.=vcov(Boot(object))} may return a full rank matrix, or it will a matrix with 
    NA corresponding to aliased regressors if same regressors are aliased in every bootstrap
    sample.
}

  \item{\dots}{arguments to pass down.}
}

\value{
Invisibly returns \code{object} for a data object, or summary for a model object.
}

\note{
The method \code{brief.matrix} calls \code{brief.data.frame}.
}

\references{

  Fox, J. and Weisberg, S. (2019)
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{
John Fox \email{jfox@mcmaster.ca}
}

\seealso{\code{\link{S}}}

\examples{
brief(rnorm(100))
brief(Duncan)
brief(OBrienKaiser, elided=TRUE)
brief(matrix(1:500, 10, 50))
brief(lm)

mod.prestige <- lm(prestige ~ education + income + type, Prestige)
brief(mod.prestige, pvalues=TRUE)
brief(mod.prestige, ~ type)
mod.mroz <- glm(lfp ~ ., data=Mroz, family=binomial)
brief(mod.mroz)
}

\keyword{manip}
