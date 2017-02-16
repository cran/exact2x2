\name{uncondExact2x2}
\alias{uncondExact2x2}
\alias{uncondExact2x2Pvals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unconditional exact tests for 2x2 tables
}
\description{
The \code{uncondExact2x2} function tests 2x2 tables assuming two independent binomial responses. Unlike the conditional exact tests which condition on both margins of the 2x2 table (see \code{\link{exact2x2}}), these unconditional tests only condition on one margin of the 2x2 table (i.e., condition on the sample sizes of the binomial responses). This makes the calculations difficult because now there is a nuisance parameter and calculations must be done over nearly the entire nuisance parameter space. 
}
\usage{
uncondExact2x2(x1, n1, x2, n2, 
    parmtype = c("difference", "ratio", "oddsratio"), nullparm = NULL, 
    alternative = c("two.sided","less", "greater"),  
    conf.int = FALSE, conf.level = 0.95, 
    method = c("FisherAdj", "simple", "score","wald-pooled", "wald-unpooled",  "user", 
      "user-fixed"), 
    tsmethod = c("central","square"), midp = FALSE, 
    gamma = 0, EplusM=FALSE, tiebreak=FALSE,
    plotprobs = FALSE, control=ucControl(), Tfunc=NULL,...)

uncondExact2x2Pvals(n1, n2, ...)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{number of events in group 1}
  \item{n1}{sample size in group 1}
  \item{x2}{number of events in group 2}
  \item{n2}{sample size in group 2}
  \item{parmtype}{
type of parameter of interest, one of "difference", "ratio" or "oddsratio" (see details)
}
  \item{nullparm}{
value of the parameter of interest at null hypothesis, NULL defaults to 0 for parmtype='difference' and 1 for parmtype='ratio' or 'oddsratio' 
}
\item{alternative}{
alternative hypothesis, one of "two.sided", "less", or "greater", default is "two.sided" (see details)
}
  \item{conf.int}{logical, calculate confidence interval?}
  \item{conf.level}{
confidence level
}
  \item{method}{
method type, one of "FisherAdj" (default), "simple", "simpleTB", "wald-pooled", "wald-unpooled", "score", "user", or "user-fixed" (see details)
}
  \item{tsmethod}{two-sided method, either "central" or "square" (see details)}
\item{midp}{logical. Use mid-p-value method?
}
\item{gamma}{Beger-Boos adjustment parameter. 0 means no adjustment. (see details).
}
\item{EplusM}{logical, do the E+M adjustment? (see details)}
\item{tiebreak}{logical, do tiebreak adjustment? (see details)}
\item{plotprobs}{logical, plot probabilities?}
\item{control}{list of algorithm parameters, see \code{\link{ucControl}}
}
\item{Tfunc}{test statistic function for ordering the sample space when method='user', ignored otherwise (see details) 
}
\item{\dots}{extra arguments passed to Tfunc (for uncondExact2x2), or passed to uncondExact2x2 (for uncondExact2x2Pvals)}
}
\details{
The \code{uncondExact2x2} function gives unconditional exact tests and confidence intervals for two independent binomial observations. The \code{uncondExact2x2Pvals} function repeatedly calls \code{uncondExact2x2} to get the p-values
for the entire sample space. 



Let X1 be binomial(n1,theta1) and X2 be binomial(n2,theta2). The parmtype determines the parameter of interest: `difference' is theta2 - theta1, 'ratio' is theta2/theta1, and `oddsratio' 
is (theta2*(1-theta1))/(theta1*(1-theta2)). 

The options \code{method}, \code{parmtype}, \code{tsmethod}, \code{alternative},
\code{EplusM}, and \code{tiebreak} define some built-in test statistic
function, Tstat, that is used to order the sample space, using \code{\link{pickTstat}} and \code{\link{calcTall}}. The first 5 arguments of Tstat must be
\code{Tstat(X1,N1,X2,N2, delta0)}, where X1 and X2 must allow vectors,
and delta0 is the null parameter value (but delta0 does not need to be used in the ordering). 
Ordering when \code{parmtype="ratio"}
or \code{parmtype="oddsratio"} is only used when there is information about the 
parameter. So the ordering function value is not used for ordering when x1=0 and x2=0 for \code{parmtype="ratio"}, and 
it is not used when (x1=0 and x2=0) or (x1=n1 and x2=n2) for \code{parmtype="oddsratio"}.


We describe the ordering functions first for the basic case,  the case when \code{tsmethod="central"} or \code{alternative!="two.sided"},
\code{EplusM=FALSE}, and \code{tiebreak=FALSE}. In this basic case the ordering function, Tstat, is determined by \code{method} and \code{parmtype}:

\itemize{
\item method='simple' - Tstat essentially replaces theta1 with x1/n1 and theta2 with x2/n2 in the parameter definition. If parmtype=`difference' 
then \code{Tstat(X1,N1,X2,N2,delta0)} returns  \code{X2/N2-X1/N1-delta0}. If parmtype='ratio' then the Tstat function
returns \code{log(X2/N2) - log(X1/N1) - log(delta0)}.
If parmtype='oddsratio' we get
\code{log( X2*(N1-X1)/(delta0*X1*(N2-X2)))}. 
\item method='wald-pooled' - Tstat is a Z statistic on the difference  using the pooled variance (not allowed if \code{parmtype!="difference"})
\item  method='wald-unpooled' - Tstat is a Z statistics on the difference using unpooled variance (not allowed if \code{parmtype!="difference"})
\item method='score' - Tstat is a Z statistic formed using score statistics,
where the parameter is defined by parmtype,
and the constrained maximum likelihood estimates of the parameter are calculated 
by \code{\link{constrMLE.difference}}, \code{\link{constrMLE.ratio}}, 
or \code{\link{constrMLE.oddsratio}}. 
\item method='FisherAdj' - Tstat is a one-sided Fisher's 'exact' mid p-value. The mid p-value is an adjustment for ties that technically removes the 'exactness' of the Fisher's p-value...BUT, here we are only using it to order the sample space, so the results of the resulting unconditional test will still be exact. 
\item method='user' - Tstat is a user supplied statistic given by \code{Tfunc}, it must be a function with the first 5 elements of its call being (X1, N1, X2, N2, delta0). The function must returns a vector of length the same as X1 and X2, where higher values suggest larger theta2 compared to theta1 (when \code{tsmethod!="square"}) or 
higher values suggest more extreme (when \code{tsmethod=="square"} and 
\code{alternative=="two.sided"}). A slower algorithm that does not require monotonicity of one-sided p-values with respect to delta0 is used.
\item method='user-fixed' - For advanced users. Tstat is a user supplied statistic given by \code{Tfunc}. It should have first 5 elements as described above but its result should not change with delta0 and it must meet Barnard's 
convexity conditions. If these conditions are met (the conditions 
are not checked, since checking them will slow the algorithm), then the p-values will be monotonic in delta0 (the null parameter for a two-sided test) and we can use a faster algorithm.
}

In the basic case, if \code{alternative="two.sided"}, the argument \code{tsmethod}="central" gives the two-sided central method. The p-value is just twice the minimum of the 
one-sided p-values (or 1 if the doubling is greater than 1). 

Now consider cases other than the basic case. 
The \code{tsmethod="square"} option gives the square of the test statistic
(when method="simple", "score", "wald-pooled", or "wald-unpooled") and larger values suggest rejection in either direction (unless method='user', then the user supplies any test statistic for which larger values suggest rejection).

The \code{tiebreak=TRUE} option breaks ties in a reasonable way when 
\code{method="simple"} (see 'details' section of \code{\link{calcTall}}).
The \code{EplusM=TRUE} option performs Lloyd's (2008) E+M ordering 
on Tstat (see 'details' section of \code{\link{calcTall}}).

If \code{tiebreak=TRUE} and \code{EplusM=TRUE}, the tiebreak calculations are always done first. 

Berger and Boos (1994) developed a very general method for calculating p-values when a nuisance parameter is present.
First, calculate a (1-gamma) confidence interval for the nuisance parameter, check for the supremum over the union of the null hypothesis parameter space 
and that confidence interval, then add back gamma to the p-value. This adjustment is valid (in other words, applied to exact tests it still gives an adjustment that is exact). The Berger-Boos adjustment is applied when \code{gamma}>0.



When method='simple' or method='user-fixed' does a simple grid search algorithm using \code{\link{unirootGrid}}.
No checks are done on the Tstat function when method='user-fixed' to make sure the simple grid search will converge to the 
proper answer. So method='user-fixed' should be used by advanced users only. 

When \code{midp=TRUE} the mid p-value is calculated (and the associated confidence interval if \code{conf.int=TRUE}) instead of the standard p-value. Loosely speaking, the standard p-value calculates the probability of observing equal or more extreme responses, while the  mid p-value calculates the probability of more extreme responses plus 1/2 the probability of equally extreme responses. The tests and confidence intervals when 
\code{midp=TRUE} are not exact, but give type I error rates and coverage of confidence intervals closer to the nominal values.
The mid p-value was studied by Lancaster (1961), see vignette on mid p-values for details. 
}
\value{
a list of class 'htest' with elements:
\item{statistic}{proportion in sample 1}
\item{parameter}{proportion in sample 2}
\item{p.value}{p-value from test}
\item{conf.int}{confidence interval on parameter given by parmtype}
\item{estimate}{MLE estimate of parameter given by parmtype}
\item{null.value}{null hypothesis value of parameter given by parmtype}
\item{alternative}{alternative hypothesis}
\item{method}{description of test}
\item{data.name}{description of data}
}
\references{
Berger, R. L. and Boos, D. D. (1994). P values maximized
over a confidence set for the nuisance parameter. Journal
of the American Statistical Association 89 1012-1016.

Lancaster, H.O. (1961). Significance tests in discrete distributions. JASA 56: 223-234.

Lloyd, C. J. (2008). Exact p-values for discrete models obtained
by estimation and maximization. Australian & New
Zealand Journal of Statistics 50 329-345.



}
\author{
Michael P. Fay, Sally A. Hunsberger
}
%\note{%%  ~~further notes~~}

\section{Warning }{
The algorithm for calculating the p-values and confidence intervals is based on a series of grid searches.
Because the grid searches are often trying to optimize non-monotonic functions, the algorithm is not 
guaranteed to give the correct answer. At the cost of increasing computation time, 
better accuracy can be obtained by increasing 
control$nPgrid,  and less often by increasing control$nCIgrid.
} 

\seealso{ See \code{\link{boschloo}} for unconditional exact tests with ordering
function based on Fisher's exact p-values. }
\examples{
# default uses method="FisherAdj"
uncondExact2x2(1,10,9,10, 
               parmtype="ratio")
uncondExact2x2(1,10,9,10, 
               method="score",parmtype="ratio")





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}