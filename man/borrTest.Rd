\name{borrTest}
\alias{borrTest}
\alias{borrPvals}
\alias{borrOrdering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Boundary-Optimized Rejection Region Test}
\description{
An unconditional exact test for the two-sample binomial problem 
when it is expected that theta1 (probability of an event in group 1)   
will be close to 1. Used for test versus control when all controls are 
expected to fail.
}
\usage{
borrTest(x1, n1, x2, n2, tuningParm = 0.025,
    parmtype = c("ratio", "difference", "oddsratio"), 
    nullparm = NULL, alternative = c("less", "greater", "two.sided"), 
    conf.int = TRUE, conf.level = 0.975,  
    controlUC = ucControl(), controlborr = borrControl(), ...)

borrPvals(n1,n2, tuningParm=0.025,    
    parmtype = c("ratio", "difference","oddsratio"), 
    nullparm = NULL, alternative = c("less", "greater","two.sided"),  
    conf.int = TRUE, conf.level = 0.975,
    controlUC=ucControl(), controlborr=borrControl(),...)

borrOrdering(n1,n2,tuningParm = .025,
            controlborr=borrControl())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
number of events in group 1
}
  \item{n1}{
sample size in group 1
}
  \item{x2}{
number of events in group 2
}
  \item{n2}{
sample size in group 2
}
  \item{tuningParm}{
tuning parameter, default is 0.025 and designs BORR tests with maximum power
for one-sided 0.025 tests
}
  \item{parmtype}{
parameter type, either 'ratio' for theta2/theta1, 
'difference' for theta2-theta1, or
'oddsratio' for theta2*(1-theta1)/(theta1*(1-theta2)).
}
  \item{nullparm}{
null parameter value, default=NULL gives parameter value for 
theta1=theta2 (e.g., 1 for 'ratio' or 0 for 'difference' ).
}
  \item{alternative}{
alternative hypothesis, BORR tests are designed for alternative='less' 
(see Note for other alternatives)
}
  \item{conf.int}{
logical, should confidence interval be calculated?
}
  \item{conf.level}{
confidence level, default is 0.975 (see note)
}
  \item{controlUC}{
a list of control parameters to define algorithms in the call to \code{\link{uncondExact2x2}}, see \code{\link{ucControl}} 
}
  \item{controlborr}{
a list of control parameters to define algorithms, see \code{\link{borrControl}}
}
  \item{\dots}{
extra arguments passed (for future development)
}
}
\details{
The boundary-optimized rejection region test is designed to test the one-sided alternative that 
theta2 < theta1, where X1 is binomial(n1,theta1), and X2 is binomial(n2,theta2). The test is 
designed to be optimal when theta1 is very close to 1. For example, in a vaccine malaria challenge study
where we expect all n1 individuals that got the control vaccine to have the event (get malaria when challenged with malaria). For details see Gabriel et al (2018). 


The function \code{borrTest} tests the results of one study, and returns 
an \code{htest} object. The function \code{borrPvals} calculates the p-values for every possible result of a study. The function \code{borrOrdering} orders every possible result of the study. 
See \code{\link{borrOrderingInternal}} for calculation details.

}
\value{
The function \code{borrPvals} returns a (n1+1) by (n2+1) matrix of p-values for all possible x1 and x2 values. The function \code{borrOrdering} returns a matrix with the rank of all possible x1 and x2 values. The function \code{borrTest} returns a list of class \code{htest} with elements: 
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
Gabriel, EE, Nason, M, Fay, MP, and Follmann, DA. (2018). A boundary-optimized rejection region test for the two-sample binomial problem. Statistics in Medicine. 37(7)  (DOI: 10.1002/sim.7579).

Antonio Martin Andres. Letter to the editor about Gabriel et al. Statistics in Medicine (to appear).

Gabriel, EE, Nason, M, Fay, MP, and Follmann, DA. Reply to letter from Martin Andres. Statistics in Medicine (to appear).
}
\author{
Martha Nason, Erin Gabriel, Michael P. Fay
}
\note{The tests are designed to have good power for the one-sided test that H0: theta2 \eqn{\ge} theta1, with 
alternative H1: theta2 < theta1 at significance level equal to \code{tuningParm}. Since the default \code{tuningParm} is 0.025, the default confidence level is 0.975 so that the confidence intervals will be compatible with the test where the one-sided p-values reject at level 0.025. 

Sometimes you may want two-sided confidence intervals on the 
parameter of interest. If you ask for a two-sided alternative, then the confidence interval and the resulting p-value will be two-sided as well. The default is a 'central' interval, so the two-sided p-value should be twice the minimum of the one-sided p-values. Further, with a \code{conf.level} of 0.95 for the two-sided alternative, the error on each side will be bounded by 0.025. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
\dontrun{borrTest(4,4,1,4)}
# Note Figure 2 in Gabriel et al is incorrect. The correct value 
# is in the response letter, and given by 
borrOrdering(4,4,tuningParm=0.025)$rankMat
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}% use one of  RShowDoc("KEYWORDS")
