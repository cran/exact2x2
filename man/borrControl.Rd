\name{borrControl}
\alias{borrControl}
\docType{data}
\title{
Algorithm variables used by borrTest.
}
\description{
Function that gives list of algorithm variables used by 
\code{\link{borrTest}}.
}
\usage{
borrControl(nAlphaGrid=10000,nThetaGrid=1000, maxIter=0, digits=4, orderFunc=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nAlphaGrid}{
number used for defining grid for searching over the (0,1) space for significance levels.
Used in \code{\link{borrOrderingAlphaGrid}}.
}
  \item{nThetaGrid}{
number of evenly spaced grid elements for searching over the (0,1) space for theta. 
Used in both \code{\link{borrOrderingAlphaGrid}} and  and \code{\link{borrOrderingByRR}}.
}
\item{maxIter}{
maximum number of searches over the alpha space. Used in  \code{\link{borrOrderingAlphaGrid}}.
}
\item{digits}{
number of digits for rounding alpha star values. Used in \code{\link{borrOrderingByRR}}.
}
\item{orderFunc}{
character vector to determine function to do the borr ordering. \code{NULL} checks for precalculated values then does \code{\link{borrOrderingAlphaGrid}}. 'AlphaGrid' uses \code{\link{borrOrderingAlphaGrid}}
and 'ByRR' uses \code{\link{borrOrderingByRR}}.
}
}
\details{
In \code{\link{borrOrderingAlphaGrid}}
we create a grid for searching over the significance level space, for the first calculation (zeroth iteration) we use \code{alpha.seq} where 
\code{alpha.seq <-
    sort(unique(c(
      10 ^ seq(log10(minAlpha), 0, length = nAlphaGrid / 2),
      seq(minAlpha, 10 ^ 0, length = nAlphaGrid / 2)
    )))}, where \code{minAlpha} is the one-sided p-value at the point (x1=n1, x2=0)
    given by \code{minAlpha<- dbinom(0,nT,nC/(nC+nT))*dbinom(nC,nC,nC/(nC+nT))}.
 If there are ties and \code{maxIter} is greater than 0, then replace each tied value with an equally spaced grid (with nAlphaGrid elements) between the adjacent non-tied values. 
 If the lowest value in the grid, \code{minAlpha}, is tied, then set \code{minAlpha<-minAlpha/10} 
 at the beginning of the iteration.
  Repeat this process up to \code{maxIter} times.
  
  
See \link{borrOrderingInternal} for more details of algorithms.
}
%\source{%%  ~~ reference to a publication or URL from which the data were obtained ~~}
%\references{%%  ~~ possibly secondary sources and usages ~~}
\examples{
borrControl(nThetaGrid=10^3)
}
\keyword{internal}
