\name{ucControl}
\alias{ucControl}
\title{
Algorithm variables used by uncondExact2x2.
}
\description{
Function that gives list of algorithm variables used by 
\code{\link{uncondExact2x2}}.
}
\usage{
ucControl(nCIgrid = 500, errbound = 0, nPgrid = 100, 
   power2 = 20, maxPgridRatio = 1 - 10^-6, 
   minPgridRatio = 10^-6)
}
\arguments{
  \item{nCIgrid}{
number of elements in the grid search for the confidence interval.
}
  \item{errbound}{
Used with large sample sizes to speed calculations, only calculate univariate binomial distribution in the middle part, exclude both tails with less than errbound/2 in each tail. When errbound=0, calculate the full distributions. 
}
  \item{nPgrid}{
number of elements to search over the null parameter space.
}
  \item{power2}{
how precise to make the grid search for the confidence interval when using the `faster' algorithm (e.g., when method='user-fixed').
}
  \item{maxPgridRatio}{
maximum binomial probability for the search over the null nuisance parameter space, when parmtype='ratio' or 'oddsratio'
}
  \item{minPgridRatio}{
maximum binomial probability for the search over the null nuisance parameter space, when parmtype='ratio' or 'oddsratio'
}
}
%\details{}
\value{
A LIST of all the named elements (see arguments to call).
}
%\references{%% }
%\author{}
%\note{}
%\seealso{}
\examples{
ucControl(errbound=10^-5)
}
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{internal}% __ONLY ONE__ keyword per line
