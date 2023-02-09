\name{power2x2}
\alias{power2x2}
\alias{ss2x2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate exact power or sample size for conditional tests for two independent binomials.
}
\description{
Power is calculated by \code{power2x2} which calls \code{\link{exact2x2}} function repeatedly. Default (strict=FALSE) does not count rejections in the wrong direction. 

Sample size is calculated by \code{ss2x2} which calls \code{power2x2} repeatedly finding the lowest sample size that has at least the nominal power, using the \code{\link{uniroot.integer}} function from the \code{ssanv} package. 
}
\usage{
power2x2(p0,p1,n0,n1=NULL,sig.level=0.05,
 alternative=c("two.sided","one.sided"),paired=FALSE,
 strict=FALSE,tsmethod=NULL,nullOddsRatio=1,
 errbound=10^-6,approx=FALSE)

ss2x2(p0,p1,power=.80,n1.over.n0=1,sig.level=0.05,
 alternative=c("two.sided","one.sided"),paired=FALSE,
 strict=FALSE,tsmethod=NULL,nullOddsRatio=1,
 errbound=10^-6,print.steps=FALSE, approx=FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p0}{ true event rate in control group }
  \item{p1}{ true event rate in treatment group}
  \item{n0}{ number of observations in control group }
 \item{n1}{ number of observations in treatment group (if NULL n1=n0)}
  \item{sig.level}{ significance level (Type I error probability) }
\item{power}{minimum power for sample size calculation} 
\item{n1.over.n0}{ratio of n1 over n0, allows for non-equal sample size allocation}
\item{alternative}{character, either "two.sided" or "one.sided", one sided tests the proper direction according to p0 and p1}
  \item{strict}{use strict interpretation of two-sided test, if TRUE counts rejections in wrong direction}
  \item{tsmethod}{two.sided method, ignored if strict=FALSE, or alternative equals 'less' or 'greater'. 
       see \code{\link{exact2x2}} for details.}
  \item{nullOddsRatio}{null odds ratio value for tests}
  \item{paired}{must be FALSE, for TRUE instead use \code{\link{powerPaired2x2}}}
\item{print.steps}{logical, print steps for calculation of sample size?}
  \item{errbound}{bound on error of calculation }
\item{approx}{give sample size or power using normal approximation only}
}
\details{
Assuming X0 ~ Binomial(n0,p0) and X1 ~ Binomial(n1,p1), 
calculates the power by repeatedly calling exact2x2 and summing probability of rejection. For speed, the function does not 
calculate the very unlikely values of X0 and X1 unless errbound=0. Power is exact, but may underestimate by at most errbound.

When strict=FALSE we do not count rejections in the wrong direction. This means that we must know the direction of the rejection, so two.sided tests are calculated as one.sided tests (in the correct direction) with level equal to sig.level/2. This is like using the tsmethod='central'. 

When \code{approx}=TRUE for \code{power2x2} use a continuity corrected normal approximation (Fleiss, 1981, p. 44). For \code{ss2x2} the calculations may be slow, so use 
\code{print.steps=TRUE} to see progress. 

}
\value{
 Both \code{power2x2} and \code{ss2x2} return an object of class 'power.htest'. A list with elements
  \item{power}{power to reject}
  \item{n0}{sample size in control group}
   \item{n1}{sample size in treatment group}
\item{p0}{ true event rate in control group }
  \item{p1}{ true event rate in treatment group}
  \item{sig.level}{ Significance level (Type I error probability) }
\item{alternative}{alternative hypothesis}
\item{note}{note about error bound}
\item{method}{description}
}
\references{ 
Fleiss. JL (1981) Statistical Methods for Rates and Proportions (second edition). Wiley. }
\author{ Michael P. Fay }
\note{ The calculations in ss2x2 can be slow when p0 is close to p1 and/or the power is large. If p0 and p1 are close with large power, it may be safer to first calculate ss2x2 with approx=TRUE
to see what the starting value will be close to. If the starting sample sizes are large (>100), it may take a while.

Note when strict=FALSE (default), the two.sided results at the 0.05 level for Fisher's exact test are like the one.sided Fisher's exact test at the 0.025 level.  
}
\section{Warning }{
There may be convergence issues using strict=FALSE with tsmethod="minlike" or "blaker" since the power is not guaranteed to be increasing in the sample size.  

} 
\seealso{See  \code{\link{ss.nonadh}} function (refinement="Fisher.exact") from the \code{ssanv} package for calculation that accounts for nonadherence in proportion of subjects. That function calls \code{\link{fisher.test}}. For power for McNemar-like test see \code{\link{powerPaired2x2}} }
\examples{
power2x2(.2,.8,12,15)
# calculate sample size with 2:1 allocation to groups
ss2x2(.2,.8,n1.over.n0=2,power=.8,approx=TRUE)
ss2x2(.2,.8,n1.over.n0=2,power=.8,print.steps=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
