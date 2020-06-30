\name{uncondPower2x2}
\alias{Power2x2}
\alias{uncondPower2x2}
\alias{SS2x2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate power or sample size for any 2x2 test.
}
\description{
The function \code{\link{Power2x2}} and \code{SS2x2} calculates the power or sample size for any 2x2 test, while the function \code{\link{uncondPower2x2}} calculates power for only tests supported by 
\code{\link{uncondExact2x2Pvals}}.
}
\usage{
Power2x2(n1, n2, theta1, theta2, alpha, pvalFunc, \dots)

uncondPower2x2(n1,n2, theta1, theta2, alpha, \dots)

SS2x2(theta1, theta2, alpha, pvalFunc, power=0.90, 
  n1start=10, increaseby=1, n2.over.n1=1,  
  maxiter=50, printSteps=TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n1}{sample size in group 1}
  \item{n2}{sample size in group 2}
  \item{theta1}{
probability of success in group 1
}
  \item{theta2}{
probability of success in group 2
}
  \item{alpha}{
significance level
}
  \item{pvalFunc}{
function that inputs x1,n1,x2,n2 and outputs a p-value.
}
\item{power}{target power}
\item{n1start}{ value of n1 for first iteration}
\item{increaseby}{positive integer, how much to increase n1 by for each iteration}
\item{n2.over.n1}{ratio of n2/n1}
\item{maxiter}{maximum number of iterations}
\item{printSteps}{logical, should the power and sample size be printed after each iteration?}
  \item{\dots}{
arguments passed to \code{\link{uncondExact2x2Pvals}} 
(for \code{uncondPower2x2}), or to \code{Power2x2} (for \code{SS2x2}). Not used and saved for future 
use for \code{Power2x2}. 
}
}
\details{
The function \code{Power2x2} is a very simple function
to calculate power. It calculates power where rejection is when the p-value from \code{pvalFunc} is less than or equal to \code{alpha}. The function \code{SS2x2} repeatedly calls \code{Power2x2} as it increases the sample size, stopping when 
the power is greater than \code{'power'}.

The function \code{\link{uncondPower2x2}} is similar except the 
p-values are calculated by \code{\link{uncondExact2x2Pvals}}.


}
\value{
the power functions return only the power. The sample size function returns a list of class 'htest.power'.
}
%\references{}
%\author{}
%\note{}

\seealso{
For power and sample size for conditional exact tests (e.g., Fisher's exact tests) see  \code{\link{power2x2}} and \code{\link{ss2x2}}. For power for the boundary-optimized rejection region (BORR) test  see \code{\link{powerBorr}}.
}
\examples{
library(exact2x2)
Power2x2(3,4,.1,.9,0.025, pvalFunc=
  function(x1,n1,x2,n2){
      boschloo(x1,n1,x2,n2, alternative="greater", 
        or=1,tsmethod="central", midp=TRUE)$p.value
  }
)
## 
\dontrun{
SS2x2(.1,.9,0.025, n1start=5, pvalFunc=
           function(x1,n1,x2,n2){
             boschloo(x1,n1,x2,n2, alternative="greater", 
                      or=1,tsmethod="central", midp=TRUE)$p.value
           }
)
}
}
\keyword{htest}% use one of  RShowDoc("KEYWORDS")