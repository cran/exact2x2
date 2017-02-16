\name{plotT}
\alias{plotT}
\alias{plotT.function}
\alias{plotT.numeric}
\alias{orderMat}
\alias{orderMat.function}
\alias{orderMat.numeric}
\docType{data}
\title{
Plot or Print  ordering function for unconditional exact test
}
\description{
The function \code{orderMat} prints the values for the ordering function for all
possible values of X1 and X2 in matrix form.

The function \code{plotT} plots the ranking of the ordering function on an n1+1 by n2+1 grid,
where each square represents a possible values for (x1,x2).
The default colors are from dark blue (highest) to light blue to white 
(middle) to light red to dard red (lowest), with black=NA.
}
\usage{
plotT(x, ...)

\method{plotT}{function}(x, n1, n2, delta0 = 1, main = "",...)


\method{plotT}{numeric}(x, n1, n2, delta0 = 1, main = "",...)


orderMat(x, ...)

\method{orderMat}{function}(x,n1,n2,delta0,graphStyle=FALSE,...)

\method{orderMat}{numeric}(x,n1,n2,delta0,graphStyle=FALSE,...)

}
\arguments{
  \item{x}{
object, either a Tstat function, or a vector of all (n1+1)*(n2+1) possible values of the function (see details).
}
  \item{n1}{sample size in group 1}
  \item{n2}{sample size in group 2}
  \item{delta0}{
null value of parameter (if needed for Tstat function)
}
  \item{main}{
plot title
}
\item{graphStyle}{
logical, order rows with lowest x1 value on the bottom?   
}
  \item{\dots}{
arguments to be passed to the Tstat function
}
}


\details{
If \code{x} is all the values of the Tstat function, then the values should be 
ordered by cycling through the x1 values (0 to n1) for each x2 value.
Specifically, it should be the result of \code{Tstat(X1,n1,X2,n2,delta0)}
where X1=rep(0:n1,n2+1) and X2=rep(0:n2,each=n1+1).

}
%\value{}
%\references{%% ~put references to the literature/web site here ~}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}
%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{

parorig<- par(no.readonly=TRUE)
par(mfrow=c(2,2),mar=c(1,3,3,1))
TT1<-pickTstat(method="score", parmtype="ratio", tsmethod="central", alternative="two.sided")
round(orderMat(TT1,8,8,1,graphStyle=TRUE),2)
TT2<-pickTstat(method="simple", parmtype="ratio", tsmethod="central", alternative="two.sided")
TT3<-pickTstat(method="simple", parmtype="difference", tsmethod="central", alternative="two.sided")
plotT(TT2, 8,8, 1, main="Ratio, Simple")
plotT(TT3, 8,8, 0, main="Difference, Simple")
plotT(TT1, 8,8, 1, main="Ratio, Score (delta0=1)")
TF<-pickTstat(method="FisherAdj", parmtype="ratio", tsmethod="central", alternative="two.sided")
plotT(TF,8,8,1, main="FisherAdj")
par(parorig)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
