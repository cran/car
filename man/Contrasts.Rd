\name{Contrasts}       
\alias{Contrasts}
\alias{contr.Treatment}
\alias{contr.Sum}
\alias{contr.Helmert}

\title{Functions to Construct Contrasts}
\description{
    These are substitutes for similarly named functions in the \pkg{stats} package 
    (note the uppercase letter starting the second word in each function name).
    The only difference is that the contrast functions from the \pkg{car} package
    produce easier-to-read names for the contrasts when they are used in statistical models. 
    
    The functions and this documentation are adapted from the \pkg{stats} package.
    }

\usage{
contr.Treatment(n, base = 1, contrasts = TRUE)

contr.Sum(n, contrasts = TRUE)

contr.Helmert(n, contrasts = TRUE)
}

\arguments{
  \item{n}{a vector of levels for a factor, or the number of levels.}
  \item{base}{an integer specifying which level is considered the baseline level. 
    Ignored if \code{contrasts} is \code{FALSE}.}
  \item{contrasts}{a logical indicating whether contrasts should be computed.}
}

\details{
    These functions are used for creating contrast matrices for use in fitting analysis of variance and regression models. 
    The columns of the resulting matrices contain contrasts which can be used for coding a factor with \code{n} levels. 
    The returned value contains the computed contrasts. If the argument \code{contrasts} is \code{FALSE} then a square matrix is returned.
    
    Several aspects of these contrast functions are controlled by options set via the \code{options} command:
    \describe{
        \item{\code{decorate.contrasts}}{This option should be set to a 2-element character vector containing the prefix and suffix
            characters to surround contrast names. If the option is not set, then \code{c("[", "]")} is used. For example, setting
            \code{options(decorate.contrasts=c(".", ""))} produces contrast names that are separated from factor names by a period.
            Setting \code{options( decorate.contrasts=c("", ""))} reproduces the behaviour of the R base contrast functions.} 
        \item{\code{decorate.contr.Treatment}}{A character string to be appended to contrast names to signify treatment contrasts;
            if the option is unset, then \code{"T."} is used.}
        \item{\code{decorate.contr.Sum}}{Similar to the above, with default \code{"S."}.}
        \item{\code{decorate.contr.Helmert}}{Similar to the above, with default \code{"H."}.} 
        \item{\code{contr.Sum.show.levels}}{Logical value: if \code{TRUE} (the default if unset), 
            then level names are used for contrasts; if \code{FALSE}, then numbers are used, as in \code{contr.sum}
            in the \code{base} package.} 
        }
        
    Note that there is no replacement for \code{contr.poly} in the \code{base} package (which produces 
    orthogonal-polynomial contrasts) since this function already constructs easy-to-read contrast names.
}

\value{
    A matrix with \code{n} rows and \code{k} columns, with \code{k = n - 1} if \code{contrasts} is \code{TRUE} 
    and \code{k = n} if \code{contrasts} is \code{FALSE}.
}

\references{
  Fox, J. and Weisberg, S. (2019) 
  \emph{An R Companion to Applied Regression}, Third Edition, Sage.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\seealso{\code{\link[stats]{contr.treatment}}, \code{\link[stats]{contr.sum}}, 
  \code{\link[stats]{contr.helmert}}, \code{\link[stats]{contr.poly}} }

\examples{
# contr.Treatment vs. contr.treatment in the base package:

lm(prestige ~ (income + education)*type, data=Prestige, 
    contrasts=list(type="contr.Treatment"))

##  Call:
##  lm(formula = prestige ~ (income + education) * type, data = Prestige,
##      contrasts = list(type = "contr.Treatment"))
##  
##  Coefficients:
##          (Intercept)                  income               education  
##              2.275753                0.003522                1.713275  
##          type[T.prof]              type[T.wc]     income:type[T.prof]  
##              15.351896              -33.536652               -0.002903  
##      income:type[T.wc]  education:type[T.prof]    education:type[T.wc]  
##              -0.002072                1.387809                4.290875  

lm(prestige ~ (income + education)*type, data=Prestige, 
    contrasts=list(type="contr.treatment"))    

##  Call:
##  lm(formula = prestige ~ (income + education) * type, data = Prestige,
##      contrasts = list(type = "contr.treatment"))
##  
##  Coefficients:
##      (Intercept)              income           education  
##          2.275753            0.003522            1.713275  
##          typeprof              typewc     income:typeprof  
##          15.351896          -33.536652           -0.002903  
##      income:typewc  education:typeprof    education:typewc  
##          -0.002072            1.387809            4.290875      
}

\keyword{models}
\keyword{regression}
