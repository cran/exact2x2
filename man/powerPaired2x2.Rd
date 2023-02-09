\name{powerPaired2x2}
\alias{powerPaired2x2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power for exact McNemar's test
}
\description{
Calculate the power for the exact McNemar's test (i.e., \code{exact2x2} with \code{paired=TRUE}) given the number of pairs
and the probability of a positive response only in the test individual in the pair (pb), and the probability of a positive response only in the control individual in the pair (pc).
}
\usage{
powerPaired2x2(pb, pc, npairs, sig.level = 0.05, 
   alternative = c("two.sided", "one.sided"), 
   strict = FALSE, nullOddsRatio = 1, errbound = 10^-6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pb}{
probability of a (0,1) response for a pair, meaning negative response in the control individual and a positive response in the test individual
}
  \item{pc}{
probability of a (1,0) response for a pair, meaning positive response in the control individual and a negative response in the test individual
}
  \item{npairs}{
the number of pairs
}
  \item{sig.level}{
significance level (also called alpha-level)
}
  \item{alternative}{
either 'one.sided' or 'two.sided' (see \code{tsmethod} for two-sided method)
}
  \item{strict}{
use strict interpretation in two-sided case (i.e., TRUE allows rejections in the 'wrong' direction)
}
  \item{nullOddsRatio}{
null odds ratio, internally passed to \code{or} argument of \code{exact2x2} with \code{paired=TRUE}
}
  \item{errbound}{
error bound, \code{errbound=0} does exact calculation, when \code{errbound>0} then speed up calculations by not calculating outcomes at either extreme with tail probabilities less than errbound/2 which may underestimate power by at most errbound.
}
\item{\dots}{
arguments passed to \code{\link{exact2x2}} (except these arguments cannot be passed this way: 
     \code{or}, \code{alternative}, \code{conf.int}, \code{paired}, \code{plot})
}
}
\details{
When \code{alternative='one.sided'} then the test automatically picks the side that is most powerful.
At this point there is no ssPaired2x2 function.
}
\value{
An object of class 'power.htest' with elements:
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{power}{power}
\item{npairs}{number of pairs}
\item{pb}{probability of a (control,test)=(0,1) response for a pair}
\item{pc}{probability of a (control,test)=(1,0) response for a pair}
\item{sig.level}{significance level or alpha-level}
\item{alternative}{either one-sided or two-sided}
\item{nullOddsRatio}{null odds ratio (or boundary between null and alternative for one-sided tests)}
\item{note}{notes about calculation (e.g., errbound value)}
\item{method}{description of method}
%% ...
}
%\references{%% ~put references to the literature/web site here ~}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
powerPaired2x2(.5,.3,npairs=20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
\keyword{htest}
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
