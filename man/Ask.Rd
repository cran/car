\name{Ask}
\alias{Ask}

\title{Change Argument to a Function Interactively}
\description{
  \code{Ask} allows you to change the argument to a function interactively.
  It is meant to be used, in lieu of a graphical control such as a slidebar,
  to adjust plotting parameters, which are
  most naturally passed as the argument to an anonymous function that sets
  up the plot.
}
\usage{
Ask(arg, fun, ...)
}

\arguments{
  \item{arg}{argument to \code{fun} to change. By specifying a vector of values,
    you can change several parameters via an argument to an anonymous function.}
  \item{fun}{function to call; often an anonymous function that sets up
    a call to plotting functions.}
  \item{\dots}{other arguments to \code{fun}; not necessary if \code{fun}
    is an anonymous function.}
}
\details{
  \code{Ask} repeatedly prompts in the R Console for the value of \code{arg}.
  To exit, enter a blank line.
}
\value{
  \code{Ask} returns invisibly the value of the last call to \code{fun};
  usually this will be \code{NULL}, and in any event is probably not of
  interest. If it is, use \code{print(Ask(arg, fun, ...))}.
}

\author{John Fox \email{jfox@mcmaster.ca}}

\examples{
    \dontrun{
attach(UN)

# enter the power-transformation parameter
# start with 1
Ask(p, function(p) qq.plot(box.cox(gdp, p), 
        ylab=paste("transformed gdp, power =",p)))

# enter an expression that evaluates to a 2-vector
# of powers; e.g., start with c(1,1); then interactively
# identify points in each plot
Ask(p, function(p) scatterplot(box.cox(gdp,p[1]), 
    box.cox(infant.mortality, p[2]), 
    xlab=paste("transformed GDP/capita, power =",p[1]),
    ylab=paste("transformed infant mortality, power =",p[2]),
    labels=rownames(UN)))    
    }
}

\keyword{iplot}
\keyword{misc}
