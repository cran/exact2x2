\name{exact2x2Plot}
\alias{exact2x2Plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot p-value function for one 2 by 2 table.
}
\description{
Plots two-sided p-values as a function of odds ratios. Can plot three types of p-values: the two-sided Fisher's exact, the central 
Fisher's exact (i.e., twice the one-sided Fisher's exact), and Blaker's exact.
}
\usage{
exact2x2Plot(x, OR = NULL, ndiv = 1000, method = "minlike", orRange = c(0.01, 10), dolog = TRUE, dolines = FALSE, dopoints = FALSE, newplot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix representing the 2 by 2 table}
  \item{OR}{odds ratio values for plot, if NULL divides orRange into ndiv pieces}
  \item{ndiv}{number of pieces to divide up odds ratio range}
  \item{method}{either "minlike","blaker" or "central"}
  \item{orRange}{range for calculating odds ratios}
  \item{dolog}{logical,plot odds ratios on log scale?}
  \item{dolines}{logical, add lines to a plot?}
  \item{dopoints}{logical, add points to a plot?}
  \item{newplot}{logical,start a new plot?}
  \item{\dots}{values passed to plot, points, or lines statement}
}
%\details{}
%\value{}
%\references{%% ~put references to the literature/web site here ~}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}
\seealso{\code{\link{exact2x2}}}
\examples{
example1<-matrix(c(6,12,12,5),2,2,dimnames=list(c("Group A","Group B"),c("Event","No Event")))
example1
exact2x2Plot(example1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
