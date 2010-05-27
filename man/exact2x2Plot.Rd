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
exact2x2Plot(x, y=NULL, OR = NULL, ndiv = 1000, tsmethod=NULL, method = NULL, paired=FALSE, 
orRange = NULL, dolog = TRUE, dolines = FALSE, dopoints = TRUE, 
doci=TRUE,  alternative=c("two.sided","less","greater"),
conf.level=.95, alphaline=TRUE, newplot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix representing the 2 by 2 table}
  \item{y}{a factor object; ignored if \code{x} is a matrix.}
  \item{OR}{odds ratio values for plot, if NULL divides orRange into ndiv pieces}
  \item{ndiv}{number of pieces to divide up odds ratio range}
  \item{tsmethod}{either "minlike","blaker" or "central"}
  \item{method}{same as tsmethod, kept for backward compatability}
\item{paired}{logical, do paired analysis giving McNemar's test p-values}
  \item{orRange}{range for calculating odds ratios}
  \item{dolog}{logical,plot odds ratios on log scale?}
  \item{dolines}{logical, add lines to a plot?}
  \item{dopoints}{logical, add points to a plot?}
  \item{doci}{logical, add vertical lines at confidence interval?}
  \item{alternative}{one of "two.sided","less","greater", type of alternative for p-values}
  \item{conf.level}{when doci=TRUE, level for confidence interval to be plotted}
  \item{alphaline}{logical, if doci=TRUE should a line be drawn at the significance level?}  
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
## add lines from central Fisher's exact
exact2x2Plot(example1,method="central",dolines=TRUE,newplot=FALSE,col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
