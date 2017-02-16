\name{constrMLE.difference}
\alias{constrMLE.difference}
\alias{constrMLE.ratio}
\alias{constrMLE.oddsratio}
\alias{constrMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate constrained MLEs.
}
\description{Calculate the constrained maximum likelihood estimate from 
2 independent binomials for the null hypothesis parameter (difference, ratio, or odds ratio of the two binomial parameters).}
\usage{
constrMLE.difference(X1, N1, X2, N2, delta0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X1}{vector, number of events in group 1}
  \item{N1}{sample size in group 1}
  \item{X2}{vector, number of events in group 2}
  \item{N2}{sample size in group 2}
  \item{delta0}{null parameter value}
}
\details{


For details see 
Farrrington and Manning (1990) for the difference, 
Miettinen and Nurminen (1985) for the ratio,
and 
Agresti and Min (2002) for the odds ratio.
}
\value{
a list with the constrained MLE parameters, p1 and p2.
}
\references{
Agresti and Min 2002, Biostatistics 3:379-386.

Farrrington and Manning, Stat in Med 1990, 1447-1454.

Miettinen, 0. and Nurminen, M. (1985). Comparative analysis
of two rates. Statistics in Medicine 4, 213-226.



}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
