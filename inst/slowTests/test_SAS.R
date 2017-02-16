

context("SAS: Barnard's test")


test_that("Barnard Exact Test p-values",{
  # 
  expect_equivalent(
    round(
      uncondExact2x2(10,35,24,61, alternative="greater", parmtype="difference", method="wald-pooled",
                     control=ucControl(nPgrid=200),conf.int=FALSE)$p.value,4),0.2142)
  # Note that we needed to go up to nPgrid>2000 to get the agreement on the last digit
  # This suggests that SAS has a better algorithm
  expect_equivalent(
    round(
      uncondExact2x2(10,35,24,61, alternative="two.sided", tsmethod="square", parmtype="difference", method="wald-pooled",
                     control=ucControl(nPgrid=5000),conf.int=FALSE)$p.value,4),0.4063)
  
  # 
  expect_equivalent(
    round(
      uncondExact2x2(4,12,9,11, alternative="greater", parmtype="difference", method="wald-pooled",
                     control=ucControl(nPgrid=100),conf.int=FALSE)$p.value,4),0.0120)
  # Note that we needed to go up to nPgrid>2000 to get the agreement on the last digit
  # This suggests that SAS has a better algorithm
  expect_equivalent(
    round(
      uncondExact2x2(4,12,9,11, alternative="two.sided", tsmethod="square", parmtype="difference", method="wald-pooled",
                     control=ucControl(nPgrid=100),conf.int=FALSE)$p.value,4),0.0215)
  
})


context("SAS: risk difference CIs")


test_that("difference CIs",{
 
  
  # Default uses central simple (inverting two one-sided using unstandardized difference)
  expect_equivalent(
    round(
      uncondExact2x2(10, 35, 24, 61, alternative="two.sided", tsmethod="central", parmtype="difference", method="simple",
                     control=ucControl(nPgrid=100),conf.int=TRUE)$conf.int,4),c(-0.1016, 0.3079))
  
  # Option: Riskdiff(method=score)
  expect_equivalent(
    round(
      uncondExact2x2(10, 35, 24, 61, alternative="two.sided", tsmethod="central",parmtype="difference", method="score",
                     control=ucControl(nPgrid=100),conf.int=TRUE)$conf.int,4),c(-0.0992, 0.2949))
  
  
  # Default uses central simple (inverting two one-sided using unstandardized difference)
  expect_equivalent(
    round(
      uncondExact2x2(4,12,9,11, alternative="two.sided", tsmethod="central", parmtype="difference", method="simple",
                     control=ucControl(nPgrid=100),conf.int=TRUE)$conf.int,4),c(0.0683, 0.7923))

  # Option: Riskdiff(method=score)
  expect_equivalent(
    round(
      uncondExact2x2(4,12,9,11, alternative="two.sided", tsmethod="central",parmtype="difference", method="score",
                     control=ucControl(nPgrid=100),conf.int=TRUE)$conf.int,4),c(0.0616, 0.7968))
  
})


context("SAS: risk ratio CIs")


test_that("ratio CIs",{
  
  
  # Default uses central simple (inverting two one-sided using unstandardized difference)
  expect_equivalent(
    round(
      uncondExact2x2(10, 35, 24, 61, alternative="two.sided", tsmethod="central", parmtype="ratio", method="simple",
                     control=ucControl(nPgrid=100),conf.int=TRUE)$conf.int,4),c(0.0409, 59024.7625))
  
  # Option: Riskdiff(method=score)
  expect_equivalent(
    round(
      uncondExact2x2(10, 35, 24, 61, alternative="two.sided", tsmethod="central",parmtype="ratio", method="score",
                     control=ucControl(nPgrid=100),conf.int=TRUE)$conf.int,4),c(0.7663, 3.0005))
  
  
  # Default uses central simple (inverting two one-sided using unstandardized difference)
  expect_equivalent(
    round(
      uncondExact2x2(4,12,9,11, alternative="two.sided", tsmethod="central", parmtype="ratio", method="simple",
                     control=ucControl(nPgrid=100),conf.int=TRUE)$conf.int,4),c(0.0795, Inf))
  
  # Option: Riskdiff(method=score)
  expect_equivalent(
    round(
      uncondExact2x2(4,12,9,11, alternative="two.sided", tsmethod="central",parmtype="ratio", method="score",
                     control=ucControl(nPgrid=100),conf.int=TRUE)$conf.int,4),c(1.0964, 10.8647))
  
})


