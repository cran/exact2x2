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
    midp=FALSE, nmc=0, eps=10^-8)
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
  \item{midp}{logical, do mid-p version of p-value and confidence intervals?}
  \item{nmc}{integer,  number of Monte Carlo replications for p-value and CI calculations, 0 (default) means calculate by numeric integration instead}
  \item{eps}{small number used to adjust numeric integration (see note)}
}
\details{
Assume X1~ Binomial(n1,p1) and X2~Binomial(n2,p2). We want to test hypotheses on a function of p1 and p2. The functions are given by parmtype: difference tests p2-p1,
ratio tests p2/p1, and odds ratio tests p2(1-p1)/(p1(1-p2)).  Let g(p1,p2) be one of the three functions. So when alternative is "less" we test H0: g(p1,p2) >= nullparm vs. H1: g(p1,p2)<nullparm.

For details when \code{midp=FALSE} see Fay, Proschan, and Brittain (2015). 

When \code{midp=TRUE}, the method performs the mid-p version on the p-value and the associated confidence intervals. 
This means that we replace the confidence distribution random variables in the p-value and CI calculations
with a random variable that is a mixture of the lower and upper CD random variables. For example, if W1L and W1U are the 
lower and upper confidence distribution random variables for group 1, then we replace those values in all calculations 
with W1midp = U1*W1L + (1-U1)*W1U, where U1 is a Bernoulli with parameter 0.5. 
For a discussion of mid-p p-values and the associated confidence intervals in a closely related context, see the vignette on mid p-values or Fay and Brittain (2016, especially the Appendix).
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
Fay, MP, Proschan, MA, and Brittain, E (2015) Combining One Sample Confidence Procedures for Inferences in the Two Sample Case. Biometrics 71: 146-156. 

Fay, Michael P., and Erica H. Brittain. (2016). Finite sample pointwise confidence intervals for a survival distribution with right-censored data. Statistics in medicine. 35: 2726-2740.
}
\author{Michael P. Fay}
\note{For numeric integration, the integrate function may have problems if nearly all of the integrand values are about 0 within the range of integration. Because of this, 
we use the eps value to make sure we integrate over ranges in which the integrand is nontrivially greater than 0. We restrict the range then add eps back to the p-value 
so that if the integrate function works perfectly, then the p-values would be very slightly conservative (for very small eps). 
There is no need to adjust the eps value. See code for detailed description of how eps is used in the calculation before changing it from the default.

An alternative method of calculation is to use Monte Carlo simulation (option with \code{nmc>0}). 
This provides a check of the numeric integration. 
There is no need to do Monte Carlo simulations for routine use. Please inform the package maintainer if the p-values or confidence intervals are substantially different when \code{nmc=0} and \code{nmc=10^7}.

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