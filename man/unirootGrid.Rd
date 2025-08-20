\name{unirootGrid}
\alias{unirootGrid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to find a root by grid search. 
}
\description{
Find the root (value where the function equals 0) of a monotonic function, \code{func}, 
using a halving algorithm grid search. 
}
\usage{
unirootGrid(func, power2 = 12, step.up = TRUE, pos.side = FALSE, 
    print.steps = FALSE, power2grid = power2gridRatio, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{func}{
monotonic function
}
  \item{power2}{
positive integer, number of grid points is 1+2^power2
}
  \item{step.up}{ logical, 
start the search at the lower end of the grid and step up?
}
  \item{pos.side}{
logical, should the root be on the positive side? In other words, should func(root)>=0
}
  \item{print.steps}{
logical, should each step that is evaluated be printed?
}
  \item{power2grid}{
  function that returns the grid. Take one argument, \code{power2}
}
  \item{\dots}{
arguments passed to \code{func}
}
}
\details{
The grid is defined with the \code{power2grid} argument that defines a function 
with an argument \code{power2}, and returns a grid with \code{1+2^power2}
elements. The root is found by a halving algorithm on the grid, so \code{func}
is calculated only \code{power2+1} times. The `root' is the element that is closest to the root, 
either on the positive side (\code{pos.side=TRUE}) or not.


The \code{unirootGrid} function calls \code{\link[ssanv]{uniroot.integer}} and finds roots based on grid search. 
The functions \code{\link[exact2x2]{power2gridRatio}} and \code{\link[exact2x2]{power2gridDifference}}
create grids for searching (0,Inf) and (-1,1) respectively. 
The \code{\link[exact2x2]{power2gridRatio}} grid is equally spaced on the log scale with about half of the grid between 0.5 and 2.
The function \code{\link[exact2x2]{power2grid}} allows more flexibility in defining grids.
}
\value{
A list with elements:

\item{iter}{number of iterations}
\item{f.root}{value of func at root}
\item{root}{root, element on the grid that is closest to the root on the negative side (if pos.side=FALSE)}
\item{bound}{interval for the  accuracy}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Michael P. Fay
}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[stats]{uniroot}} and \code{\link[ssanv]{uniroot.integer}}
}
\examples{
# print.steps prints all iterations, 
# with x=rank of grid value (e.g., x=1 is lowest value in grid) 
# f(x) really is f(grid[x]) where grid is from the power2grid function 
unirootGrid(function(x){ x - .37}, power2=10, power2grid=power2gridRatio, 
  print.steps=TRUE, pos.side=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{optimize}% use one of  RShowDoc("KEYWORDS")
%\keyword{internal}% __ONLY ONE__ keyword per line
