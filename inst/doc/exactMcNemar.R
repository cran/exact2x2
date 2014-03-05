### R code from vignette source 'exactMcNemar.Rnw'

###################################################
### code chunk number 1: exactMcNemar.Rnw:23-24
###################################################
library(exact2x2)


###################################################
### code chunk number 2: exactMcNemar.Rnw:87-89
###################################################
x<-matrix(c(21,9,2,12),2,2)
mcnemar.test(x)


###################################################
### code chunk number 3: exactMcNemar.Rnw:92-93
###################################################
mcnemar.test(x,correct=FALSE)


###################################################
### code chunk number 4: exactMcNemar.Rnw:123-124
###################################################
mcnemar.exact(x)


