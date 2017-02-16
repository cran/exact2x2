\name{pickTstat}
\alias{pickTstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pick T statistic (ordering function) for unconditional exact test.}
\description{
Called from \code{\link{uncondExact2x2}}.
}
\usage{
pickTstat(method, parmtype, tsmethod, alternative)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
method type, one of "simple", "wald-pooled", "wald-unpooled", "score" (see details) 
}
\item{parmtype}{
type of parameter of interest, one of "difference", "ratio" or "oddsratio" 
(see details)}
\item{tsmethod}{two-sided method, either "central" or "square" (see details) }
  \item{alternative}{
alternative hypothesis, one of "two.sided", "less", or "greater", default is "two.sided" 
}
}
\details{
See 'details' section of \code{\link{uncondExact2x2}}.
}
\value{
A function that has the following arguments:
  \item{X1}{vector, number of events in group 1}
  \item{N1}{sample size in group 1}
  \item{X2}{vector, number of events in group 2}
  \item{N2}{sample size in group 2}
  \item{delta0}{null parameter value}
  
and outputs a vector the same length as X1 and X2.  
}
%\references{%% ~put references to the literature/web site here ~}
%\author{}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
