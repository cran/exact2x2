\name{exact2x2-package}
\alias{exact2x2-package}
\docType{package}
\title{
Exact Tests and Confidence Intervals for 2x2 Tables
}
\description{
There are 7 main functions in the package. 
The \code{\link{exact2x2}} function calculates the exact conditional tests with matching confidence intervals as detailed in Fay (2010a <DOI:10.1093/biostatistics/kxp050>,2010b). The functions \code{\link{ss2x2}} and \code{\link{power2x2}} calculate the sample size and power related to the tests of \code{exact2x2}. The \code{\link{uncondExact2x2}}
and \code{\link{boschloo}} functions calculate unconditional exact tests.
The \code{\link{binomMeld.test}} function calculates melded confidence intervals for two sample binomial inferences (see Fay, Proschan, and Brittain, 2015 <DOI:10.1111/biom.12231>).
Finally, the \code{\link{borrTest}} function calculates the boundary optimized rejection region test that 
creates unconditional exact tests that have power optimized when group 1 is expected to have 100 percent failure. For example, in vaccine challenge studies where the control group are all expected to get infected (see Gabriel, et al, 2018 <DOI:10.1002/sim.7579>, the letter about that paper by Martin Andres <DOI:10.1002/sim.7630>, and the response <DOI:10.1002/sim.7684>).
 }
\details{
\tabular{ll}{
Package: \tab bpcp\cr
Type: \tab Package\cr
Version: \tab 1.6.4.1\cr
Date: \tab 2020-06-29\cr
License: \tab GPL3\cr
LazyLoad: \tab yes\cr
}
}
\author{
Michael P. Fay, Sally A. Hunsberger, Martha Nason, Erin Gabriel

Maintainer: Michael P. Fay <mfay@niaid.nih.gov>

}
\references{
Fay, M. P. (2010a). Confidence intervals that Match Fisher's exact and Blaker's exact tests. Biostatistics,
  11: 373-374 (go to doc directory for earlier version or 
  \url{https://www.niaid.nih.gov/about/brb-staff-fay} for link to official version).

Fay, M.P. (2010b). Two-sided Exact Tests and Matching Confidence Intervals for Discrete Data. R Journal 2(1):53-58.


Fay, MP, Proschan, MA, and Brittain, E (2015). Combining One Sample Confidence Procedures for Inference in the Two Sample Case. Biometrics. 71: 146-156. 

Gabriel, EE, Nason, M, Fay, MP, and Follmann, DA. (2018). A boundary-optimized rejection region test for the two-sample binomial problem. Statistics in Medicine. 37(7): 1047-1058  (DOI: 10.1002/sim.7579).


Gabriel, EE, Nason, M, Fay, MP, and Follmann, DA. (2018). Reply to letter from Martin Andres. Statistics in Medicine 37(14): 2303-2306.

Martin Andres, Antonio. (2018). Letter to the editor about Gabriel et al. Statistics in Medicine 37(14) 2301-2302.

}
\keyword{ package }
\keyword{htest}
\keyword{nonparametric}
%\seealso{ }
%\examples{}
