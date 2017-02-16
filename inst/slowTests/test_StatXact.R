library(exact2x2)
library(testthat)

context("Testing StatXact version 11: Difference, Invert 2 1-sided, gamma=0")


test_that("central",{
out<-uncondExact2x2(2,7,6, 13,parmtype="difference", tsmethod="central", method="score", conf.int=TRUE)
# get p-values and conf intervals from StatXact 11 and check
expect_equal(round(c(out$p.value,out$conf.int),4),c(0.5569,-0.3108,0.5787))
out<-uncondExact2x2(0,7,13, 13,parmtype="difference", tsmethod="central", method="score", conf.int=TRUE)
# The CIs below agree with SAS Ver 9.4, SAS calculates p-value as Barnard's it is p<0.0001  
expect_equal(round(c(out$p.value,out$conf.int),c(9,4,4)),c(4.758e-06,0.5904,1.000))
out<-uncondExact2x2(2,7,11, 13,parmtype="ratio", tsmethod="central", method="score", conf.int=TRUE)
expect_equal(round(c(out$p.value,out$conf.int),c(5,3,1)),c(0.02064,1.078,29.1))
})


context("Testing StatXact version 11: Invert 2-sided, gamma=0")


test_that("square",{
  out<-uncondExact2x2(2,7,6, 13,parmtype="difference", tsmethod="square", method="score", conf.int=TRUE)
  # get p-values and conf intervals from StatXact 11 and check
  expect_equal(round(c(out$p.value,out$conf.int),c(4,4,4)),c(0.5374,-0.3079,0.5542))
  out<-uncondExact2x2(0,7,13, 13,parmtype="difference", tsmethod="square", method="score", conf.int=TRUE)
  # Following answers are from StatXact Version 11, fails because StatXact appears to be wrong
  expect_equal(round(c(out$p.value,out$conf.int),c(8,4,4)),c(2.44e-06,0.6088,0.6101))
  out<-uncondExact2x2(2,7,11, 13,parmtype="ratio", tsmethod="square", 
                      method="score", conf.int=TRUE)
  ## does not equal StatXact 11 CI: 
  expect_equal(round(c(out$p.value,out$conf.int),c(4,2,2)),c(0.0134,1.14,14.45))
})

context("Testing StatXact version 1:  Invert 2 1-sided, gamma=1e-06")


test_that("central",{
  out<-uncondExact2x2(2,7,6, 13,parmtype="difference", tsmethod="central", method="score", gamma=1e-06, conf.int=TRUE)
  # get p-values and conf intervals from StatXact 11 (answers do not change from gamma=0) and check
  expect_equal(round(c(out$p.value,out$conf.int),4),c(0.5569,-0.3108,0.5787))
  out<-uncondExact2x2(0,7,13, 13,parmtype="difference", tsmethod="central", method="score", gamma=1e-06, conf.int=TRUE)
  # The CIs below agree with SAS Ver 9.4, SAS calculates p-value as Barnard's it is p<0.0001  
  expect_equal(round(c(out$p.value,out$conf.int),c(9,4,4)),c(6.758e-06,0.5904,1.000))
  out<-uncondExact2x2(2,7,11, 13,parmtype="ratio", tsmethod="central", method="score", gamma=1e-06, conf.int=TRUE)
  expect_equal(round(c(out$p.value,out$conf.int),c(5,3,2)),c(0.02064,1.078,28.74))
})

context("Testing StatXact version 11: difference only with gamma=0 different from gamma=1e-06")


test_that("central example 2:",{
  pgam<-uncondExact2x2(10,35, 24, 61,parmtype="difference", tsmethod="central", method="score", gamma=1e-06, conf.int=TRUE) 
  # get p-values and conf intervals from StatXact 11 
  expect_equal(round(c(pgam$p.value,pgam$conf.int),c(4,5,4)),c(0.3012,-0.09731, 0.2949))
  p<-uncondExact2x2(10,35, 24, 61,parmtype="difference", tsmethod="central", method="score", gamma=0, conf.int=TRUE)$p.value 
  # get p-values and conf intervals from StatXact 11 and check
  expect_equal(round(c(pgam$p.value,pgam$conf.int),c(4,5,4)),c(0.4284,-0.09922, 0.2949))
})


context("Testing StatXact version 11: Boschloo")

test_that("central example 2:",{
  pgr<-boschloo(10,35, 24, 61, alternative="greater", conf.int=FALSE)$p.value 
  # get p-values  from StatXact 11 
  expect_equal(round(c(pgr),c(4)),c(0.1513))
  pml<-boschloo(10,35, 24, 61, alternative="two.sided", tsmethod="minlike", conf.int=FALSE)$p.value 
  # get p-value from StatXact 11 and check
  expect_equal(round(c(pml),c(4)),c(0.3381))
})

