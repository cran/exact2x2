\name{exact2x2}
\alias{exact2x2}
\alias{fisher.exact}
\alias{blaker.exact}
\alias{mcnemar.exact}
%\alias{exact2x2-package}
\title{Exact Conditional Tests for 2 by 2 Tables of Count Data}
\description{
  Performs exact conditional tests for two by two tables. For independent binary responses, performs 
either Fisher's exact test or Blaker's exact test for testing hypotheses about the odds ratio. 
The commands follow the style of  \code{\link[stats]{fisher.test}}, the difference is that
for two-sided tests there are three methods for calculating the exact test, and for each of the three methods
its matching 
confidence interval is returned (see details).  
For paired binary data resulting in a two by two table, performs an exact McNemar's test.}
\usage{
exact2x2(x, y = NULL, or = 1, alternative = "two.sided",
    tsmethod = NULL, conf.int = TRUE, conf.level = 0.95,
    tol = 0.00001, conditional = TRUE, paired=FALSE, 
    plot=FALSE, midp=FALSE)
fisher.exact(x, y = NULL, or = 1, alternative = "two.sided",
    tsmethod = "minlike", conf.int = TRUE, conf.level = 0.95,
    tol = 0.00001, midp=FALSE)
blaker.exact(x, y = NULL, or = 1, alternative = "two.sided",
    conf.int = TRUE, conf.level = 0.95, tol = 0.00001)
mcnemar.exact(x,y=NULL, conf.level=.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{either a two-dimensional contingency table in matrix form,
    or a factor object.}
  \item{y}{a factor object; ignored if \code{x} is a matrix.}
  \item{or}{the hypothesized odds ratio. Must be a single numeric.  }
  \item{alternative}{indicates the alternative hypothesis and must be
    one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
    if "two.sided" uses method defined by tsmethod.}
  \item{tsmethod}{one of "minlike","central", or "blaker". NULL defaults to "minlike" when 
  paired=FALSE and "central" when paired=TRUE or midp=TRUE.
Defines type of two-sided method (see details). Ignored if alternative="less" or "greater".}
  \item{conf.int}{logical indicating if a confidence interval should be
    computed.}
  \item{conf.level}{confidence level for the returned confidence
    interval.  Only used if
    \code{conf.int = TRUE}.}
  \item{tol}{tolerance for confidence interval estimation.}
  \item{conditional}{TRUE. Unconditional exact tests should use \code{\link{uncondExact2x2}}.}
  \item{paired}{logical. TRUE gives exact McNemar's test, FALSE are all other tests}
\item{midp}{logical. TRUE gives mid p-values and mid-p CIs. Not supported for tsmethod='minlike' or 'blaker' }
\item{plot}{logical. TRUE gives basic plot of point null odds ratios by p-values, for greater plot control use \code{\link{exact2x2Plot}}. Not supported for midp=TRUE.}
}
\details{
  The motivation for this package is to match the different two-sided conditional exact tests 
  for 2x2 tables with the appropriate confidence intervals. 

  There are three ways to calculate the two-sided conditional exact tests, 
  motivated by three different ways to define the p-value. 
  The usual two-sided Fisher's exact test defines the p-value as the sum of probability
  of tables with 
  smaller likelihood than the observed table (\code{tsmethod}="minlike"). 
The central Fisher's exact test defines the p-value as twice the one-sided p-values
(but with a maximum p-value of 1). Blaker's (2000) exact test defines the p-value
as the sum of the tail probibility in the observed tail plus the largest tail probability
in the opposite tail that is not greater than the observed tail probability. 

   In \code{\link[stats]{fisher.test}} the p-value uses the two-sample method 
associated with \code{tsmethod}="minlike", but the confidence interval method 
associated with  \code{tsmethod}="central". The probability that the
 lower central confidence limit is less than the true odds ratio is bounded by 
\code{1-(1-conf.level)/2} for the central intervals, but not for the other two two-sided 
methods.   
 The confidence intervals in for \code{exact2x2} match the test associated 
with alternative. In other words, the confidence interval is the smallest interval that contains the confidence set that is 
the inversion of the associated test (see Fay, 2010).
The functions \code{fisher.exact} and \code{blaker.exact} are just wrappers for certain
options in  \code{exact2x2}.

  If \code{x} is a matrix, it is taken as a two-dimensional contingency
  table, and hence its entries should be nonnegative integers.
  Otherwise, both \code{x} and \code{y} must be vectors of the same
  length.  Incomplete cases are removed, the vectors are coerced into
  factor objects, and the contingency table is computed from these.

  P-values are obtained directly using the (central or non-central) hypergeometric
  distribution. 

  The null of conditional
  independence is equivalent to the hypothesis that the odds ratio
  equals one.  \sQuote{Exact} inference can be based on observing that in
  general, given all marginal totals fixed, the first element of the
  contingency table has a non-central hypergeometric distribution with
  non-centrality parameter given by the odds ratio (Fisher, 1935).  The
  alternative for a one-sided test is based on the odds ratio, so
  \code{alternative = "greater"} is a test of the odds ratio being bigger
  than \code{or}.

 When paired=TRUE, this denotes there is some pairing of the data. For example, 
 instead of Group A and Group B, we may have pretest and posttest binary responses. 
 The proper two-sided test for such a setup is McNemar's Test, which only uses the off-diagonal
 elements of the 2x2 table, and tests that both are equal or not. The exact version 
 is based on the binomial distribution on one of the off-diagonal values conditioned on the total 
 of both off-diagonal values. We use \code{\link[exactci]{binom.exact}} from the \code{exactci} package, and convert the 
 p estimates and confidence intervals (see note) to odds ratios (see Breslow and Day, 1980, p. 165). The function 
 \code{mcnemar.exact} is just a wrapper to call \code{exact2x2} with \code{paired=TRUE, alternative="two.sided",tsmethod="central"}. 
 One-sided exact McNemar-type tests may be calculated using the  \code{exact2x2} function with \code{paired=TRUE}.
For details of McNemar-type tests see Fay (2010, R Journal). 

The mid p-value is an adjusted p-value to account for discreteness. The mid-p adjustment is not guaranteed to give type I error rates that are less than or equal to nominal levels, but gives p-values that lead to the probability of rejection that is sometimes less than the nominal level and sometimes greater than the nominal level. This adjustment is sometimes used because exact p-values for discrete data cannot give actual type I error rates equal to the nominal value unless randomization is done (and that is not typically done because two researchers doing the same method could get different answers). Essentially, exact p-values lead to the probability of rejecting being less than the nominal level for most parameter values in the null hypothesis in order to make sure that it is not greater than the nominal level for ANY parameter values in the null hypothesis. The mid p-value was studied by Lancaster (1961), and for the 2x2 case by Hirji et al (1991).
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{p.value}{the p-value of the test}
  \item{conf.int}{a confidence interval for the odds ratio}
  \item{estimate}{an estimate of the odds ratio.  Note that the
    \emph{conditional} Maximum Likelihood Estimate (MLE) rather than the
    unconditional MLE (the sample odds ratio) is used.}
  \item{null.value}{the odds ratio under the null, \code{or}.}
  \item{alternative}{a character string describing the alternative
    hypothesis}
  \item{method}{a character string, changes depending on alternative and tsmethod}
  \item{data.name}{a character string giving the names of the data}
}
\references{
Blaker, H. (2000) Confidence curves and improved exact confidence intervals for discrete distributions.
Canadian Journal of Statistics 28: 783-798.

Breslow, NE and Day NE (1980). Staistical Methods in Cancer Research: Vol 1-The analysis of Case-Control Studies.
    IARC Scientific Publications. IARC, Lyon.

Fay, M. P. (2010). Confidence intervals that Match Fisher's exact and Blaker's exact tests. Biostatistics,
  11: 373-374 (go to doc directory for earlier version or 
  \url{https://www.niaid.nih.gov/about/brb-staff-fay} for link to official version).

Fay M.P. (2010). Two-sided Exact Tests and Matching Confidence Intervals for Discrete Data. R Journal 2(1):53-58.

Fisher, R.A. (1935) The logic of inductive inference. Journal of the Royal Statistical Society
Series A 98:39-54.

Hirji, K.F., Tan, S-J, and Elashoff, R.M. (1991). A quasi-exact test for comparing two binomial proportions. Statistics in Medicine 10: 1137-1153.

Lancaster, H.O. (1961). Significance tests in discrete distributions. JASA 56: 223-234.
}
\author{
Michael Fay
}
\note{ The default exact confidence intervals for the odds ratio when paired=TRUE (those matching the exact McNemar's test) 
are transformations of the Clopper-Pearson exact confidence intervals for a single binomial parameter which are central intervals. 
See note for \code{\link[exactci]{binom.exact}} for discussion of exact binomial confidence intervals. }
\seealso{
  \code{\link[stats]{fisher.test}} or \code{\link[stats]{mcnemar.test}}
}
\examples{
## In example 1, notice how fisher.test rejects the null at the 5 percent level, 
## but the 95 percent confidence interval on the odds ratio contains 1 
## The intervals do not match the p-value.
## In fisher.exact you get p-values and the matching confidence intervals 
example1<-matrix(c(6,12,12,5),2,2,dimnames=list(c("Group A","Group B"),c("Event","No Event")))
example1
fisher.test(example1)
fisher.exact(example1,tsmethod="minlike")
fisher.exact(example1,tsmethod="central")
blaker.exact(example1)
## In example 2, this same thing happens, for
## tsmethod="minlike"... this cannot be avoided because 
## of the holes in the confidence set.
##  
example2<-matrix(c(7,255,30,464),2,2,dimnames=list(c("Group A","Group B"),c("Event","No Event")))
example2
fisher.test(example2)
exact2x2(example2,tsmethod="minlike")
## you can never get a test-CI inconsistency when tsmethod="central"
exact2x2(example2,tsmethod="central")
}
\keyword{htest}
