\name{mcnemarExactDP}
\alias{mcnemarExactDP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exact McNemar (Paired Binary) Test with Difference in Proportions}
\description{
Gives a valid (i.e., exact) test of paired binary responses, with compatible confidence intervals on the difference in proportions.
}
\usage{
mcnemarExactDP(x, m, n, nullparm = 0, alternative = c("two.sided", "less", "greater"), 
    conf.level = 0.95, nmc = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
number of pairs with mismatched responses  
}
  \item{x}{
number of pairs with response of 1 for treatment and 0 for control
}
\item{n}{
total number of pairs
}
  \item{nullparm}{
null parameter value for the difference in proportions: proportion with events on treatment minus proportion with events on control
}
  \item{alternative}{
alternative hypothesis, must be one of "two.sided", "greater" or "less"}
  \item{conf.level}{
confidence level for the returned confidence interval
}
  \item{nmc}{
number of Monte Carlo replications, nmc=0 (default) uses numeric integration instead
}
}
\details{
For paired binary responses, a simple test is McNemars test, which conditions on the number of discordant pairs. 
The \code{\link{mcnemar.exact}} function gives results in terms of odds ratios. This function gives results in terms 
of the difference in proportions. The p-values will be identical between the two functions, but the estimates and 
confidence intervals will be different.

For this function, we use the melding idea (Fay, et al, 2015), to create compatable confidence intervals
    with exact versions of McNemars test. For details see Fay and Lumbard (2021).
 See Fagerland, et al (2013) for other parameters and methods related to paired binary responses. 
    The advantage of this version is that it is exact, and faster than the unconditional exact methods (which may be more powerful). 
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{p.value}{the p-value of the test}
  \item{conf.int}{a confidence interval for the difference in proportions}
  \item{estimate}{sample proportions and their difference}
  \item{null.value}{difference in proportions under the null}
  \item{alternative}{a character string describing the alternative
    hypothesis}
  \item{method}{a character string describing the test}
  \item{data.name}{a character string giving the names of the data}
}
\references{
Fay, MP, Proschan, MA, and Brittain, E (2015). Combining one-sample confidence procedures for inference in the two-sample case. Biometrics,71(1),146-156.

Fay MP, and Lumbard, K (2021). Confidence Intervals for Difference in Proportions for Matched Pairs Compatible with Exact McNemars or Sign Tests. Statistics in Medicine, 40(5): 1147-1159.

Fagerland, Lydersen and Laake (2013), Recommended tests and confidence
intervals for paired binomial proportions. Statitics in Medicine, 33:2850-2875.
}
\author{
Michael P. Fay, Keith Lumbard
}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{mcnemar.exact}} or \code{\link{exact2x2}} with paired=TRUE for confidence intervals on the odds ratio.
}
\examples{
# For test on contingency table of the pairs
# From Bentur, et al (2009) Pediatric Pulmonology 44:845-850.
# see also Table II of Fagerland, Lydersen and Laake 
# (2013, Stat in Med, 33: 2850-2875)
# 
#                 After SCT
#                 AHR       No AHR
#              ----------------- 
# Before SCT  |    
#         AHR |    1         1 
#      No AHR |    7        12
#              -----------------

ahr<-matrix(c(1,7,1,12),2,2,
            dimnames=list(paste("Before SCT,",c("AHR","No AHR")),
                          paste("After SCT,",c("AHR","No AHR"))))
mcnemarExactDP(n=sum(ahr),m=ahr[1,2]+ahr[2,1], x=ahr[1,2])
# compare to mcnemar.exact
# same p-value, but mcnemar.exact gives conf int on odds ratio 
mcnemar.exact(ahr)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
