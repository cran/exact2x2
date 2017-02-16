## Test the E+M option
## Test is proposed in Lloyd (2008, Aust. N.Z. J. Stat 50(4), 329-345)
## This script reproduces some results from Table 1
### the T rows (top is one-sided, bottom is two-sided)
## and the E+M column: The p-value results are 0.02518 (one-sided)  and 0.03681 (two-sided) 
library(exact2x2)

x1<-14
n1<-47
x2<-48
n2<-283


context("E+M from Lloyd 2008: Check results from Table 1")

test_that("E+M wald-pooled",{
## one-sided
  expect_equal(
      round(uncondExact2x2(x1,n1,x2,n2,parmtype = "difference", alternative="less", method="wald-pooled", 
                           EplusM=TRUE, control=ucControl(nPgrid=1000),conf.int=FALSE)$p.value,5), 0.02518)

## two-sided
  expect_equal(
      round(uncondExact2x2(x1,n1,x2,n2,parmtype = "difference", tsmethod="square", method="wald-pooled", 
                           EplusM=TRUE, control=ucControl(nPgrid=1000),conf.int=FALSE)$p.value,5),0.03681)
  
})

