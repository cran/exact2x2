\name{borrOrderingInternal}
\alias{borrOrderingInternal}
\alias{borrOrderingAlphaGrid}
\alias{borrOrderingByRR}
\alias{borrOrderingPreCalc}
\alias{borrPreCalc}
\alias{calcRejectProb}
\alias{getThreshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BORR Ordering, internal calculation functions
}
\description{
Three functions for calculating the BORR ordering. The default did some slow \code{borrOrdering} calculations are done ahead of time and stored in sysdata.rda, for n1 and n2 smaller than 21 and for \code{tuningParm=0.025}. 
}
\usage{
borrOrderingAlphaGrid(n1, n2, tuningParm = 0.025, controlborr = borrControl())

borrOrderingByRR(n1, n2, tuningParm = 0.025, controlborr = borrControl())

borrOrderingPreCalc(n1, n2, tuningParm=0.025, orderPreCalc=orderPreCalc)

borrPreCalc(NList=seq(2,20),
           tuningParm = 0.025,
           controlborr = borrControl())
           
calcRejectProb(p.ctrl, Threshold, p.trt = p.ctrl, n.trt, n.ctrl, max.uninf.ctrls = n.ctrl)

getThreshold(n.ctrl, n.trt, tuningParm = 0.025, nThetaGrid = 1000, 
    max.uninf.ctrls = n.ctrl, forceConvex = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n1}{
sample size in group 1
}
  \item{n2}{
sample size in group 2
}
  \item{tuningParm}{
tuning parameter, default is 0.025 and designs BORR tests with maximum power
for one-sided 0.025 tests
}
  \item{controlborr}{
a list of control parameters to define algorithms, see \code{\link{borrControl}}
}
\item{orderPreCalc}{
a list of precalculated orderings (see details)
}
\item{NList}{
list of n1 and n2 values for creating \code{orderPreCalc} object. Does all possible combinations
}
\item{p.ctrl}{
vector of theta values for theta1, usually determined by \code{controlborr$nThetaGrid} 
}
\item{Threshold}{
vector of threshold values that define one rejection region.
}
\item{p.trt}{
vector of theta values for theta2, usually determined by \code{controlborr$nThetaGrid}
}
\item{n.trt}{
n2 (notation matches the Gabriel, et al paper)
}
\item{n.ctrl}{
n1 (notation matches the Gabriel, et al paper)
}
\item{max.uninf.ctrls}{
set to n.ctrl, see code before changing it
}
\item{forceConvex}{
logical, should always be TRUE. If you want to try FALSE check the code first.
}
}
\details{
All BORR ordering functions automatically enforce Barnard's convexity in the 
rejection regions (in response to the letter of Martin Andres). 
Note that the original ordering in Figure 2 of Gabriel et al was incorrect. 
The correct value is in the response  letter by Gabriel et al (see also the example code 
in \code{\link{borrTest}}).


The \code{controlborr$orderFunc} determines which function calculates the borr ordering. 
When  \code{controlborr$orderFunc=NULL} (the default) the code first searches to see if there is a precalculated ordering (see below), and if not it calls \code{borrOrderingByRR} if \code{n1+n2<=16},
and otherwise calls \code{borrOrderingAlphaGrid}. 
When  \code{controlborr$orderFunc='AlphaGrid'} then it calls \code{borrOrderingAlphaGrid},
when \code{controlborr$orderFunc='ByRR'} then it calls \code{borrOrderingByRR}.


The function \code{borrOrderingByRR} calculates the ordering based on trying convex 
rejection regions and calculating the alpha star value when different points that are added are 
just barely rejected. This leads to fast and accurate calculates for small n1 and n2 (less than 8), but can be slow for larger n1 and n2. It rounds the alpha star values to the nearest \code{controlborr$digits}, to avoid computer problems with ties (remember the alpha star values themselves are calculated by a grid on the theta values).


The function \code{borrOrderingAlphaGrid} calculates the ordering based on a grid of alpha values.
It can be faster for larger n1 and n2, but its accuracy depends on the \code{controlborr$nAlphaGrid}. 

The function \code{borrPreCalc} as run in the example should produce \code{orderPreCalc}. It was actually run on a parallel processing machine as 361 separate jobs. These calculations can 
take a bit of time. Then \code{\link{borrOrdering}} (when \code{controlborr$orderFunc=NULL})
will automatically check to see if the ordering has previously been calculated and if so will call 
\code{borrOrderingPreCalc} and if not call \code{borrOrderingAlphaGrid}.

The functions \code{calcRejectProb} and \code{getThreshold} are called by both \code{borrOrderingAlphaGrid} and \code{borrOrderingByRR}. 
}
\value{
The function \code{borrOrderingAlphaGrid} and \code{borrOrderingByRR}  returns an rank matrix
as well as an alpha matrix. 
The alpha matrix is the minimum alpha for each point to just enter the rejection region
(in the notation of Gabriel et al, it is Min(alphastar: delta(alphastar, NC, NT, YC, YT)=1)).  
The rank matrix is the ordering matrix as in Figure 2 (see correction in letter).
The \code{borrOrderingPreCalc} only returns the rank matrix. 
The list \code{orderPreCalc} has elements:
\item{orderList}{ a list of the rank matrices, with orderList[[i]] associated with n1List[i] and n2List[i]}
\item{controlborr}{control used in calculating orderings, see \code{\link{borrControl}}}
\item{tuningParm}{the tuning parmeter used in the orderings}
\item{n1List}{the n1List used in the orderings}
\item{n2List}{the n2List used in the orderings}
}
\references{
Gabriel, EE, Nason, M, Fay, MP, and Follmann, DA. (2018). A boundary-optimized rejection region test for the two-sample binomial problem. Statistics in Medicine. 37(7)  (DOI: 10.1002/sim.7579).

Antonio Martin Andres. Letter to the editor about Gabriel et al. Statistics in Medicine (to appear).

Gabriel, EE, Nason, M, Fay, MP, and Follmann, DA. Reply to letter from Martin Andres. Statistics in Medicine (to appear).
}
%\note{}
\examples{
\dontrun{
# This is the call that should produce the orderPreCalc object 
# used by borrOrderingPreCalc
orderPreCalc<-borrPreCalc(NList=2:20,
    tuningParm = 0.025,
    controlborr = borrControl(nAlphaGrid = 10000, 
    nThetaGrid=1000, maxIter=0)) 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
