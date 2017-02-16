


library(exact2x2)

context("Odds ratio with score examples")


test_that("",{
  # central, see Agresti and Min 2002, p. 381
  out<- uncondExact2x2(1,26,2,26, alternative="two.sided", tsmethod="central", parmtype="oddsratio", method="score", 
                       control=ucControl(nPgrid=1000), conf.int=TRUE)
  expect_equal(
    round(c(out$conf.int[1],out$conf.int[2]),c(2,1)),c(0.10, 127.2)
    )
  # square, see Agresti and Min 2002, p. 381
  outs<- uncondExact2x2(1,26,2,26, alternative="two.sided", tsmethod="square", parmtype="oddsratio", method="score", 
                       control=ucControl(nPgrid=1000), conf.int=TRUE)
  expect_equal(
    round(c(outs$conf.int[1],outs$conf.int[2]),c(2,1)),c(0.15, 62.7)
  )


  # square, see Fagerland, et al 2015, SMMR, 224-254, Table 6
  outf<- uncondExact2x2(1,34,7,34, alternative="two.sided", tsmethod="square", parmtype="oddsratio", method="score", 
                        control=ucControl(nPgrid=1000), conf.int=TRUE)
  
  
  
  
  expect_equal(
    round(c(outf$conf.int[1],outf$conf.int[2]),c(2,0)),c(1.19, 72)
  )  

})


context("Problem with Delta Close to 1: oddsratio score")


test_that("Sally found error Feb 15 2017:",{
  out<-uncondExact2x2(2,12,7, 12,
                      alternative=c("two.sided"),
                      nullparm=NULL,
                      parmtype=c("oddsratio"),
                      conf.int=TRUE,
                      conf.level=.95,
                      method=c("score"),
                      tsmethod=c("central")  )
  # get p-values and conf intervals from StatXact 11 and check
  expect_equal(round(c(out$p.value,out$conf.int),c(3,2,0)),c(0.045, 1.040, 103))
})


