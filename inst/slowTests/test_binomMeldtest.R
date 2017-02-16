# test binomMeld.test

library(exact2x2)


context("binomMeld.test: midp=FALSE difference")


test_that("",{
  set.seed(121)
  bmc<-binomMeld.test(20,20,20,20, parmtype="difference", nullparm=0, alternative="greater", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="difference", nullparm=0, alternative="greater", midp=FALSE, nmc=0)
  summaryBM<-function(b,rdigits=c(7,7,7)){
      c(round(b$p.value,rdigits[1]), round(b$conf.int[1],rdigits[2]), round(b$conf.int[2],rdigits[3]))    
  }
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)

  bmc<-binomMeld.test(20,20,20,20, parmtype="difference", nullparm=0, alternative="less", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="difference", nullparm=0, alternative="less", midp=FALSE, nmc=0)

  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  bmc<-binomMeld.test(0,20,0,20, parmtype="difference", nullparm=0, alternative="greater", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="difference", nullparm=0, alternative="greater", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  bmc<-binomMeld.test(0,20,0,20, parmtype="difference", nullparm=0, alternative="less", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="difference", nullparm=0, alternative="less", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
  set.seed(513)
  bmc<-binomMeld.test(4,20,7,12, parmtype="difference", nullparm=0, alternative="greater", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="difference", nullparm=0, alternative="greater", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(121)
  bmc<-binomMeld.test(4,20,7,12, parmtype="difference", nullparm=0, alternative="less", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="difference", nullparm=0, alternative="less", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
})


context("binomMeld.test: midp=FALSE  ratio")


test_that("",{
  set.seed(121)
  bmc<-binomMeld.test(20,20,20,20, parmtype="ratio", nullparm=1, alternative="greater", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="ratio", nullparm=1, alternative="greater", midp=FALSE, nmc=0)
  summaryBM<-function(b,rdigits=c(7,7,7)){
    c(round(b$p.value,rdigits[1]), round(b$conf.int[1],rdigits[2]), round(b$conf.int[2],rdigits[3]))    
  }
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(2)
  bmc<-binomMeld.test(20,20,20,20, parmtype="ratio", nullparm=1, alternative="less", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="ratio", nullparm=1, alternative="less", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(1)
  bmc<-binomMeld.test(0,20,0,20, parmtype="ratio", nullparm=1, alternative="greater", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="ratio", nullparm=1, alternative="greater", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  bmc<-binomMeld.test(0,20,0,20, parmtype="ratio", nullparm=1, alternative="less", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="ratio", nullparm=1, alternative="less", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
 
  
  
  set.seed(513)
  bmc<-binomMeld.test(4,20,7,12, parmtype="ratio", nullparm=1, alternative="greater", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="ratio", nullparm=1, alternative="greater", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(90)
  bmc<-binomMeld.test(4,20,7,12, parmtype="ratio", nullparm=1, alternative="less", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="ratio", nullparm=1, alternative="less", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
   
})


context("binomMeld.test: midp=FALSE  oddsratio")


test_that("",{
  set.seed(121)
  bmc<-binomMeld.test(20,20,20,20, parmtype="oddsratio", nullparm=1, alternative="greater", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="oddsratio", nullparm=1, alternative="greater", midp=FALSE, nmc=0)
  summaryBM<-function(b,rdigits=c(7,7,7)){
    c(round(b$p.value,rdigits[1]), round(b$conf.int[1],rdigits[2]), round(b$conf.int[2],rdigits[3]))    
  }
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(21)
  bmc<-binomMeld.test(20,20,20,20, parmtype="oddsratio", nullparm=1, alternative="less", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="oddsratio", nullparm=1, alternative="less", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(513)
  bmc<-binomMeld.test(0,20,0,20, parmtype="oddsratio", nullparm=1, alternative="greater", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="oddsratio", nullparm=1, alternative="greater", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(121)
  bmc<-binomMeld.test(0,20,0,20, parmtype="oddsratio", nullparm=1, alternative="less", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="oddsratio", nullparm=1, alternative="less", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
  
  set.seed(513)
  bmc<-binomMeld.test(4,20,7,12, parmtype="oddsratio", nullparm=1, alternative="greater", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="oddsratio", nullparm=1, alternative="greater", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(31)
  bmc<-binomMeld.test(4,20,7,12, parmtype="oddsratio", nullparm=1, alternative="less", midp=FALSE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="oddsratio", nullparm=1, alternative="less", midp=FALSE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
  
})

context("binomMeld.test: midp=TRUE difference")


test_that("",{
  set.seed(121)
  bmc<-binomMeld.test(20,20,20,20, parmtype="difference", nullparm=0, alternative="greater", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="difference", nullparm=0, alternative="greater", midp=TRUE, nmc=0)
  summaryBM<-function(b,rdigits=c(7,7,7)){
      c(round(b$p.value,rdigits[1]), round(b$conf.int[1],rdigits[2]), round(b$conf.int[2],rdigits[3]))    
  }
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)

  bmc<-binomMeld.test(20,20,20,20, parmtype="difference", nullparm=0, alternative="less", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="difference", nullparm=0, alternative="less", midp=TRUE, nmc=0)

  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  bmc<-binomMeld.test(0,20,0,20, parmtype="difference", nullparm=0, alternative="greater", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="difference", nullparm=0, alternative="greater", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  bmc<-binomMeld.test(0,20,0,20, parmtype="difference", nullparm=0, alternative="less", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="difference", nullparm=0, alternative="less", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
  set.seed(513)
  bmc<-binomMeld.test(4,20,7,12, parmtype="difference", nullparm=0, alternative="greater", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="difference", nullparm=0, alternative="greater", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(121)
  bmc<-binomMeld.test(4,20,7,12, parmtype="difference", nullparm=0, alternative="less", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="difference", nullparm=0, alternative="less", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
})


context("binomMeld.test: midp=TRUE  ratio")


test_that("",{
  set.seed(121)
  bmc<-binomMeld.test(20,20,20,20, parmtype="ratio", nullparm=1, alternative="greater", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="ratio", nullparm=1, alternative="greater", midp=TRUE, nmc=0)
  summaryBM<-function(b,rdigits=c(7,7,7)){
    c(round(b$p.value,rdigits[1]), round(b$conf.int[1],rdigits[2]), round(b$conf.int[2],rdigits[3]))    
  }
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(2)
  bmc<-binomMeld.test(20,20,20,20, parmtype="ratio", nullparm=1, alternative="less", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="ratio", nullparm=1, alternative="less", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(1)
  bmc<-binomMeld.test(0,20,0,20, parmtype="ratio", nullparm=1, alternative="greater", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="ratio", nullparm=1, alternative="greater", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  bmc<-binomMeld.test(0,20,0,20, parmtype="ratio", nullparm=1, alternative="less", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="ratio", nullparm=1, alternative="less", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
 
  
  
  set.seed(513)
  bmc<-binomMeld.test(4,20,7,12, parmtype="ratio", nullparm=1, alternative="greater", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="ratio", nullparm=1, alternative="greater", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(90)
  bmc<-binomMeld.test(4,20,7,12, parmtype="ratio", nullparm=1, alternative="less", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="ratio", nullparm=1, alternative="less", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
   
})


context("binomMeld.test: midp=TRUE  oddsratio")


test_that("",{
  set.seed(121)
  bmc<-binomMeld.test(20,20,20,20, parmtype="oddsratio", nullparm=1, alternative="greater", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="oddsratio", nullparm=1, alternative="greater", midp=TRUE, nmc=0)
  summaryBM<-function(b,rdigits=c(7,7,7)){
    c(round(b$p.value,rdigits[1]), round(b$conf.int[1],rdigits[2]), round(b$conf.int[2],rdigits[3]))    
  }
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(21)
  bmc<-binomMeld.test(20,20,20,20, parmtype="oddsratio", nullparm=1, alternative="less", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(20,20,20,20, parmtype="oddsratio", nullparm=1, alternative="less", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(513)
  bmc<-binomMeld.test(0,20,0,20, parmtype="oddsratio", nullparm=1, alternative="greater", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="oddsratio", nullparm=1, alternative="greater", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(121)
  bmc<-binomMeld.test(0,20,0,20, parmtype="oddsratio", nullparm=1, alternative="less", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(0,20,0,20, parmtype="oddsratio", nullparm=1, alternative="less", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
  
  set.seed(513)
  bmc<-binomMeld.test(4,20,7,12, parmtype="oddsratio", nullparm=1, alternative="greater", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="oddsratio", nullparm=1, alternative="greater", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  set.seed(31)
  bmc<-binomMeld.test(4,20,7,12, parmtype="oddsratio", nullparm=1, alternative="less", midp=TRUE, nmc=10^6)
  bint<-binomMeld.test(4,20,7,12, parmtype="oddsratio", nullparm=1, alternative="less", midp=TRUE, nmc=0)
  
  expect_equal(summaryBM(bmc), summaryBM(bint), check.attributes=FALSE, tolerance=.005)
  
  
  
})
