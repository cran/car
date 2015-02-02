\name{Soils}
\alias{Soils}
\docType{data}
\title{Soil Compositions of Physical and Chemical Characteristics}
\description{
 Soil characteristics were measured on samples from three types of
 contours (Top, Slope, and Depression) and at four depths (0-10cm,
10-30cm, 30-60cm, and 60-90cm).  The area was divided into 4 
blocks, in a randomized block design. (Suggested by Michael Friendly.)
}
\usage{Soils}
\format{
  A data frame with 48 observations on the following 14 variables.  There are 3 factors and 9 response variables.
  \describe{
    \item{\code{Group}}{a factor with 12 levels, corresponding to the combinations of \code{Contour} and \code{Depth} }
    \item{\code{Contour}}{a factor with 3 levels: \code{Depression} \code{Slope} \code{Top}}
    \item{\code{Depth}}{a factor with 4 levels: \code{0-10} \code{10-30} \code{30-60} \code{60-90}}
    \item{\code{Gp}}{a factor with 12 levels, giving abbreviations for the groups: 
	\code{D0} \code{D1} \code{D3} \code{D6} \code{S0} \code{S1} \code{S3} \code{S6} \code{T0} \code{T1} \code{T3} \code{T6}}
    \item{\code{Block}}{a factor with levels \code{1} \code{2} \code{3} \code{4}}
    \item{\code{pH}}{soil pH}
    \item{\code{N}}{total nitrogen in \%}
    \item{\code{Dens}}{bulk density in gm/cm$^3$ }
    \item{\code{P}}{total phosphorous in ppm}
    \item{\code{Ca}}{calcium in me/100 gm.}
    \item{\code{Mg}}{magnesium in me/100 gm.}
    \item{\code{K}}{phosphorous in me/100 gm.}
    \item{\code{Na}}{sodium in me/100 gm.}
    \item{\code{Conduc}}{conductivity}
  }
}
\details{
These data provide good examples of MANOVA and canonical discriminant analysis in a somewhat
complex multivariate setting.  They may be treated as a one-way design (ignoring \code{Block}),
by using either \code{Group} or \code{Gp} as the factor, or a two-way randomized block
design using \code{Block}, \code{Contour} and \code{Depth} (quantitative, so orthogonal
polynomial contrasts are useful).
}
\source{
  Horton, I. F.,Russell, J. S., and Moore, A. W. (1968)
  Multivariate-covariance and canonical analysis: 
  A method for selecting the most effective discriminators in a multivariate situation.
  \emph{Biometrics} \bold{24}, 845--858.
  Originally from \samp{http://www.stat.lsu.edu/faculty/moser/exst7037/soils.sas} but no longer available there.
}
\references{
  Khattree, R., and Naik, D. N. (2000)
  \emph{Multivariate Data Reduction and Discrimination with SAS Software.}
   SAS Institute.
  
  Friendly, M. (2006)
  Data ellipses, HE plots and reduced-rank displays for
  multivariate linear models: SAS software and examples.
  \emph{Journal of Statistical Software}, 17(6),
  \url{http://www.jstatsoft.org/v17/i06}.  

}

\keyword{datasets}
