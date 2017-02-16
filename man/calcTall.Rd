\name{calcTall}
\alias{calcTall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate all Tstat for all values of the (n1+1) X (n2+1) sample space from the two sample binomial problem.
}
\description{
Used mostly by internal call from \code{\link{uncondExact2x2}}.
If \code{EplusM=FALSE} and \code{tiebreak=FALSE} then the result is just
\code{Tstat(allx,n1,ally,n2,delta0)}. Otherwise does more complicated 
calculations.
}
\usage{
calcTall(Tstat, allx, n1, ally, n2, delta0 = 0, parmtype = "difference", 
    alternative = "two.sided", tsmethod = "central", EplusM = FALSE, tiebreak = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tstat}{
ordering function
}
  \item{allx}{
  vector of x1 values, typically rep(0:n1,n2+1)
}
  \item{n1}{
sample size in group 1 
}
  \item{ally}{
  vector of x2 values, typically rep(0:n2,each=n1+1)
}
  \item{n2}{
sample size in group 2
}
  \item{delta0}{
null parameter value for input into Tstat
}
  \item{parmtype}{
parmeter type, either 'difference', 'ratio', or 'oddsratio'
}
  \item{alternative}{
alternative hypothesis, either 'two.sided' or not
}
  \item{tsmethod}{
two-sided method, either 'central' or 'square'  
}
  \item{EplusM}{
logical, do E+M ordering of Lloyd (2008)?
}
  \item{tiebreak}{
  logical, do tie break method? Only allowed when tsmethod!='square'.
}
}
\details{


When \code{tiebreak=TRUE} does a method that breaks ties in the ordering function differently  depending on the \code{parmtype} value. The tie breaks are developed to make sense when \code{method="simple"} and \code{tsmethod!="square"}, when applied to other methods it may not necessarily break ties reasonably. For that reason \code{tiebreak=TRUE} returns an error when \code{tsmethod="square"}. For \code{parmtype="difference"} ties are broken based on Z scores on the difference in proportions, with larger values of \code{Z}
treated as larger. This means that when the sample proportions are equal,
the ties are not broken. For \code{parmtype="ratio"} ties are broken based
on \code{abs(Z)}, where the Z scores are based on the difference in log proportions, except when x1=0 (when ties are broken by x2) or x2=0 (when ties are broken by 1/x1). For \code{parmtype="oddsratio"} ties are broken based
on \code{abs(Z)}, where here the Z scores are based on the 
difference in log odds, except when x1=0 or x1=n1 or x2=0 or x2=n2 (see code for specifics). 

The E+M method, is to take an existing ordering function, Tstat, and calculate
a one-sided p-value based on that ordering function evaluated at the constrained 
maximum likelihood estimates of the parameters. The ordering is then 
the set of one-sided p-values from Pr[T(X)<=T(xobs)],  except when 
\code{alternative="two.sided"} and \code{tsmethod="square"} in which case 
it is 1-p, where p, the p-value, is based on Pr[T(X)>=T(xobs)]. The latter exception is needed so that larger values are more likely to reject.

If \code{tiebreak=TRUE} and \code{EplusM=TRUE}, the teibreak calculations are always done first. 
}
\value{
a vector of the same length as allx, giving values of Tstat function at all
values in the sample space.
}
%\references{}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
