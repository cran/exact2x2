\name{exact2x2}
\alias{exact2x2}
\alias{fisher.exact}
\alias{blaker.exact}
\alias{exact2x2-package}
\title{Exact Conditional Tests for 2 by 2 Tables of Count Data}
\description{
  Performs either Fisher's exact test or Blaker's exact test for testing hypotheses about the odds ratio from a two by two contingency 
table with fixed marginals. The commands follow the style of  \code{\link[stats]{fisher.test}}, the difference is that the 
confidence intervals are matching ones (see details).}
\usage{
exact2x2(x, y = NULL, or = 1, alternative = "two.sided", tsmethod = "minlike", conf.int = TRUE, conf.level = 0.95, tol = 0.00001, conditional = TRUE)
fisher.exact(x, y = NULL, or = 1, alternative = "two.sided", tsmethod = "minlike", conf.int = TRUE, conf.level = 0.95, tol = 0.00001)
blaker.exact(x, y = NULL, or = 1, alternative = "two.sided", conf.int = TRUE, conf.level = 0.95, tol = 0.00001)
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
  \item{tsmethod}{one of "minlike","central", or "blaker". 
Defines type of two-sided method (see details). Ignored if alternative="less" or "greater".}
  \item{conf.int}{logical indicating if a confidence interval should be
    computed.}
  \item{conf.level}{confidence level for the returned confidence
    interval.  Only used if
    \code{conf.int = TRUE}.}
  \item{tol}{tolerance for confidence interval estimation.}
  \item{conditional}{TRUE. Unconditional exact tests not supported at this time.}
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
the inversion of the associated test (see Fay, 2009).
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

Fay, M. P. (2009). Confidence intervals for Fisher's exact and Blaker's exact tests.
  (unpublished manuscript, see pdf in doc directory of this package).

Fisher, R.A. (1935) The logic of inductive inference. Journal of the Royal Statistical Society
Series A 98:39-54.
}
\author{
Michael Fay
}
%\note{}
\seealso{
  \code{\link[stats]{fisher.test}}
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
## In example 2, this same thing happens, in some examples... this cannot be avoided because of the 
## holes in the confidence set.
##  
example2<-matrix(c(7,255,30,464),2,2,dimnames=list(c("Group A","Group B"),c("Event","No Event")))
example2
fisher.test(example2)
exact2x2(example2,tsmethod="minlike")
exact2x2(example2,tsmethod="central")
exact2x2(example2,tsmethod="blaker")
}
\keyword{htest}