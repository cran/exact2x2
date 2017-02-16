\name{power2grid}
\alias{power2grid}
\alias{power2grid}
\alias{power2gridRatio}
\alias{power2gridDifference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create grid for root search.
}
\description{
Used with \code{\link{unirootGrid}}
}
\usage{
power2grid(power2 = 3, from = 10, to = 1, dolog = TRUE)
power2gridRatio(power2 = 3)
power2gridDifference(power2 = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power2}{
positive integer, determines length of grid, \code{length(grid)=1+2^power2}
}
  \item{from}{
lowest value of grid
}
  \item{to}{
highest value of grid
}
  \item{dolog}{
  logical, make grid equally spaced on the log scale?
  }
}
\details{
The function \code{power2gridRatio} gives a grid for searching from 0 to Inf
equally spaced on the log scale,
with about half of the observations from 0.5 to 2.
and \code{power2gridDifference} gives an equally spaced grid for searching from -1 to 1.
}
\value{a vector for grid search of length \code{1+2^power2} for use in
\code{\link{unirootGrid}}
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{unirootGrid}}
}
\examples{
power2gridRatio(3)
power2gridDifference(3)
power2grid(3,from=-1,to=1,dolog=FALSE)
power2grid(3, from=1,to=9, dolog=FALSE)
power2grid(3, from=1,to=9, dolog=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
