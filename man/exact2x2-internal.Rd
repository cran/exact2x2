\name{exact2x2-internal}
\alias{exact2x2CI}
\alias{mcnemar.exact.calc}
\alias{fisherCalcMidp}
%\name{fisher.blaker.pvals}
\alias{exact2x2Pvals}
\alias{binomMeldCalcInt}
\alias{binomMeldCalcMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal functions for exact2x2. Not to be called by user.}
\description{
The function \code{exact2x2Pvals} can calculate p-values for a vector of odds ratios.
The function \code{exact2x2CI} is the code that calculates the confidence intervals for the two-sided 
Fisher's exact test and Blaker's exact test.  The functions \code{binomMeldCalcInt} and \code{binomMeldCalcMC}
are called by \code{\link{binomMeld.test}}.
}
\usage{
fisherCalcMidp(x,or,alternative,conf.int,conf.level)
exact2x2Pvals(x, or, relErr=1+10^(-7),tsmethod = "minlike", alternative="two.sided")
exact2x2CI(x, tsmethod="minlike", conf.level=0.95, tol=0.00001, orRange=c(10^-10,10^10))
mcnemar.exact.calc(bb,cc,or,alternative,tsmethod="central",conf.level=.95, midp=FALSE)
binomMeldCalcInt(x1,n1,x2,n2,nullparm=NULL, 
           parmtype=c("difference","oddsratio","ratio"),
           conf.level=0.95, conf.int=TRUE,
           alternative=c("two.sided","less","greater"), midp=FALSE, eps=10^-8)
binomMeldCalcMC(x1,n1,x2,n2,nullparm=NULL, 
                           parmtype=c("difference","oddsratio","ratio"),
                           conf.level=0.95, conf.int=TRUE,
                           alternative=c("two.sided","less","greater"),
                           midp=FALSE,nmc=10^6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix representing 2 by 2 table}
  \item{or}{odds ratio, may be a vector}
  \item{relErr}{relative error. This is used to handle true ties on the computer. (see details). }
  \item{tsmethod}{either "minlike","blaker", "central"}
  \item{conf.int}{logical, calculate CI?}
  \item{conf.level}{confidence level}
  \item{tol}{tolerance}
  \item{orRange}{range for search for odds ratio confidence interval}
 \item{alternative}{indicates the alternative hypothesis and must be
    one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.}
\item{midp}{logical. Do midp adjustment?
}
}
\details{P-values for both the two-sided Fisher's exact and Blaker's exact test add probabilities from the opposite tail if 
either the cumulative probabilities (for Blaker's test) or the probabilities (Fisher's test) are less than or equal to those 
of the observed tail. Since sometimes the p-value at odds ratio=1 is important, we may have problems if the opposite tail is 
some very small different probability due to computer rounding, when mathematically the probabilities are exactly the same 
and should be included. To get around this problem \code{\link[stats]{fisher.test}} uses relErr so that it chooses all
d<= d[i]*relErr and if mathematically d[i] is equal to another value in d but there is a slightly computer rounding error, that value
will be included. We use the same tactic here.


The function \code{mcnemar.exact.calc} is just a simple call to \code{binom.test} with 
\code{p=.5}.


 }
\value{
  Output from \code{exact2x2Pvals} is a LIST, with
  \item{or}{vector of odds ratios}
  \item{pvals}{vector of two-sided p-values}
 Output from \code{exact2x2CI} is a confidence interval with attributes: conf.level and conf.int.prec (a list 
of the bounds on the precision of the limits).
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ ~Make other sections like Warning with \section{Warning }{....} ~}
\seealso{\code{\link{exact2x2}} }
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
