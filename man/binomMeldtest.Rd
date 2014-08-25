\name{binomMeld.test}
\alias{binomMeld.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Melded Binomial Confidence Intervals and Tests 
}
\description{
Creates tests to compare two binomials, giving confidence intervals for either the difference in proportions, the rate ratio, or the odds ratio. The 95 percent confidence intervals have been shown to guarantee nominal coverage by extensive numerical calculations. It has been theoretically proven that the p-values from the one-sided tests on the null hypothesis of equality match Fisher's exact p-values.  
}
\usage{
binomMeld.test(x1, n1, x2, n2, nullparm = NULL, 
    parmtype = c("difference", "oddsratio", "ratio"), 
    conf.level = 0.95, conf.int=TRUE, 
    alternative = c("two.sided", "less", "greater"), 
    eps = 10^-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{number of events in group 1}
  \item{n1}{sample size in group 1}
  \item{x2}{number of events in group 2}
  \item{n2}{sample size in group 2}
  \item{nullparm}{value of the parameter of interest at null, default of NULL gives 0 for parmtype='difference' and 1 for parmtype='ratio' or 'oddsratio'}
  \item{parmtype}{type of parameter of interest, one of "difference", "ratio" or "oddsratio" (see details)}
  \item{conf.level}{confidence level}
  \item{conf.int}{logical, calculate confidence intervals?}
  \item{alternative}{alternative hypothesis, one of "two.sided", "less", or "greater" (see details)}
  \item{eps}{small number used to adjust numeric integration (see note)}
}
\details{
Assume X1~ Binomial(n1,p1) and X2~Binomial(n2,p2). We want to test hypotheses on a function of p1 and p2. The functions are given by parmtype: difference tests p2-p1,
ratio tests p2/p1, and odds ratio tests p2(1-p1)/(p1(1-p2)).  Let g(p1,p2) be one of the three functions. So when alternative is "less" we test H0: g(p1,p2) >= nullparm vs. H1: g(p1,p2)<nullparm.

For details see Fay, Proschan, and Brittain (2014). 
}
\value{
An object of class 'htest'. A list with elements
\item{statistic }{proportion of events in group 1}
\item{parameter}{proportion of events in group 2}
\item{p.value}{p-value}
\item{conf.int}{confidence interval}
\item{estimate}{estimate of g(p1,p2) by plugging in sample proportions, i.e., unconditional MLE}
\item{null.value}{value of g(p1,p2) under null}
\item{alternative}{type of alternative hypothesis}
\item{method}{description of test}
\item{data.name}{character explicit description of data}
}
\references{
Fay, MP, Proschan, MA, and Brittain, E (2014) Combining One Sample Confidence Procedures for Inferences in the Two Sample Case (to appear in Biometrics). 

}
\author{Michael P. Fay}
\note{For numeric integration, the integrate function may have problems if nearly all of the integrand values are about 0 within the range of integration. Because of this, 
we use the eps value to make sure we integrate over ranges in which the integrand is nontrivially greater than 0. We restrict the range then add eps back to the p-value 
so that if the integrate function works perfectly, then the p-values would be very slightly conservative (for very small eps). 
There is no need to adjust the eps value. See code for detailed description of how eps is used in the calculation before changing it from the default.
}
%\section{Warning }{....} ~
%\seealso{}
\examples{
# Note the p-value for all tests of equality 
# (Null Hypthesis: true prop 1=true prop 2)
# are the same, and equal to the 
# Fisher's exact (central) p-value
binomMeld.test(3,5,1,8,parmtype="difference")
binomMeld.test(3,5,1,8,parmtype="ratio")
# note that binomMeld.test gives the unconditional MLE 
# for the odds ratio, while fisher.test and exact2x2 
# gives the conditional MLE for the odds ratio
# (also fisher.test gives the odds ratio defined as 
#  the inverse of how it is defined in binomMeld.test)
binomMeld.test(3,5,1,8,parmtype="oddsratio")
exact2x2(matrix(c(1,8-1,3,5-3),2,2),tsmethod="central")

}
\keyword{htest }