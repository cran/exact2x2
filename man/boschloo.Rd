\name{boschloo}
\alias{boschloo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boschloo's test for 2x2 Tables}
\description{
Boschloo's test is an exact unconditional test for 2x2 tables based on ordering the sample space by Fisher's exact p-values. This function generalizes that test in several ways (see details).
}
\usage{
boschloo(x1, n1, x2, n2, alternative = c("two.sided", "less", "greater"), 
    or = NULL, conf.int = FALSE, conf.level = 0.95, midp = FALSE, 
    tsmethod = c("central", "minlike"), control=ucControl())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{number of events in group 1}
  \item{n1}{sample size in group 1}
  \item{x2}{number of events in group 2}
  \item{n2}{sample size in group 2}
    \item{alternative}{
alternative hypothesis, one of "two.sided", "less", or "greater", default is "two.sided" (see details)
}
\item{or}{odds ratio under the null hypothesis}
  \item{conf.int}{logical, calculate confidence interval?}
  \item{conf.level}{
confidence level
}
\item{midp}{logical. Use mid-p-value method?
}
  \item{tsmethod}{two-sided method, either "central" or "minlike" (see details)}
\item{control}{list of algorithm parameters, see \code{\link{ucControl}}
}
}
\details{
The traditional Boschloo (1970) test is to use Fisher's exact p-values (under the null that p1=p2) to order the sample space and to use that ordering to perform an unconditional exact test. Here we generalize this to test for different null hypothesis values (other that odds ratios of 1). 

For the two-sided alternatives, the traditional method uses 
tsmethod='minlike' (for example, in the Exact R package)
but our default is tsmethod='central'. The one-sided tests use ordering by 
the appropriate p-value (or 1 minus the p-value for alternative='greater' so that the ordering function follows our convention for user supplied ordering functions, see method='user' option in \code{\link{uncondExact2x2}}).

The option \code{midp} orders the sample space by the mid-p value associated with 
Fisher's exact test, and additionally gives mid-p values. This means that unlike the \code{midp=FALSE} case, when \code{midp=TRUE} the test is not exact (i.e., guaranteed to bound the type I error rate at the nominal level), but has type I error rates that are on average (over the possible null parameter values) closer to the nominal level. 

If you want to order by the mid-p values from Fisher's exact test but get an exact test, use the \code{method="FisherAdj"} with \code{midp=FALSE} in \code{\link{uncondExact2x2}}. 

The \code{boschloo} function only gives confidence intervals for the odds ratio, for getting confidence intervals closely related to Boschloo p-values (but not exactly matching Boschloo p-values) for the difference or ratio, use \code{\link{uncondExact2x2}} with \code{method="FisherAdj"}.
}
\value{
a list of class 'htest' with elements:
\item{statistic}{proportion in sample 1}
\item{parameter}{proportion in sample 2}
\item{p.value}{p-value from test}
\item{conf.int}{confidence interval on odds ratio}
\item{estimate}{odds ratio estimate}
\item{null.value}{null hypothesis value of odds ratio}
\item{alternative}{alternative hypothesis}
\item{method}{description of test}
\item{data.name}{description of data}
}
\references{
Boschloo, R. D. "Raised conditional level of significance for the 2x2-table when testing the equality of two probabilities." Statistica Neerlandica 24.1 (1970): 1-9.
}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{   \code{exact.test} in package \CRANpkg{Exact} for Boschloo test p-value computation. Also see \code{method"FisherAdj"} in \code{\link{uncondExact2x2}} for a closely related test. }
\examples{
# defaults to the central two-sided version
boschloo(1,5,6,7)
boschloo(1,5,6,7,alternative="greater")
## traditional two-sided Boschloo test (not central!)
boschloo(1,5,6,7, tsmethod="minlike")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")