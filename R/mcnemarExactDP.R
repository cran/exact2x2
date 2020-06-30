mcnemarExactDP <- function(x,
                                  m,
                                  n,
                                  nullparm=0, 
                                  alternative=c("two.sided","less", "greater"),
                                  conf.level = 0.95,
                                  nmc = 0) {
  alternative<- match.arg(alternative)
  # check input values...
  if (length(x)!=1 | length(n)!=1 | length(m)!=1){ stop("n,m, and x must be scalars (vectors of length 1)")}
  if (x>m){ stop("x>m, not possible since x is the number of positives out of m")}
  if (m>n){ stop("m>n, not possible since m is the number of mis-matched responses out of n pairs")}
  # one-sided p-value for use in deciding whether to 
  # use TL or TU
  alpha<- 1- conf.level
  if (alternative=="two.sided"){ 
    alpha<- alpha/2 
    pU0<- pbeta(0.5,x+1,m-x, lower.tail=FALSE)
    pL0<- pbeta(0.5,x, m-x+1)
  } else if (alternative=="less"){
    pU0<- pbeta(0.5,x+1,m-x, lower.tail=FALSE)
    pL0<- NULL
  } else if (alternative=="greater"){
    pL0<- pbeta(0.5,x, m-x+1)
    pU0<-NULL
  }
  # if nmc=0 then use numeric integration
  # otherwise use Monte Carlo methods
  
  # define extremes for one-sided intervals, then fill in the other side
  lower<- -1
  upper<- 1
  
  if (nmc > 0) {
    # Monte Carlo Method
    TL <- rbeta(nmc, m, n - m + 1)
    TU <- rbeta(nmc, m + 1, n - m)
    BL <- rbeta(nmc, x, m - x + 1)
    BU <- rbeta(nmc, x + 1, m - x)
    d <- function(T, B) {
      T * (2 * B - 1)
    }
    if (alternative=="two.sided" | alternative=="greater"){
      lower <- min(quantile(d(TL, BL), probs = alpha),
                   quantile(d(TU, BL), probs = alpha))
      DL<-d(TL,BL)
      DU<-d(TU,BL)
      if (nullparm==0){
        pL<-pL0
      } else {
        pL<- max( length(DL[DL<=nullparm])/nmc, 
                  length(DU[DU<=nullparm])/nmc)
      }
    }
    if (alternative=="two.sided" | alternative=="less"){
      upper <- max(quantile(d(TL, BU), probs = 1 - alpha),
                   quantile(d(TU, BU), probs = 1 - alpha))
      DL<-d(TL,BU)
      DU<-d(TU,BU)
      if (nullparm==0){
        pU<-pU0
      } else {
        pU<- max( length(DL[DL>=nullparm])/nmc, 
                  length(DU[DU>=nullparm])/nmc)             
      }
    }
    
  } else {
    # Numeric Integration method
    pval.1sided <- function(Delta,X,M,N, alt="less",...) {
      x<-X
      m<-M
      n<-N
      if (alt=="less"){
        if (Delta==0){
          pval<- pbeta(0.5, x+1, m-x, lower.tail=FALSE)
        } else if (Delta>0){
          if (m==n){
            # m=n -> TU~Beta(m+1,0)=1, so 
            # Pr[TU(2BU-1)>=Delta] = Pr[2BU-1 >=Delta]=Pr[BU>= 1/2 + Delta/2]
            pval<-pbeta(0.5 + Delta /2, x+1, m-x, lower.tail=FALSE)
          } else if (m==x){
            # x=m => BU ~ Beta(x+1,0)=1, so Pr[TU(2BU-1)>=Delta] = Pr[TU>=Delta]
            pval<-pbeta(Delta,m+1,n-m, lower.tail=FALSE)
          } else {
            ifunc <- function(tt, delta = Delta) {
              pbeta(0.5 + delta / (2 * tt), x+1, m-x, lower.tail=FALSE) * dbeta(tt, m+1, n-m)
            }
            # do not integrate of the ifunc, unless there are at least 2 non-zero values (>1e-16)
            fvals <- ifunc(seq(0, 1, length.out = 10000))
            if(sum(fvals > 10e-16) > 1){
            nonzero <- range(which(fvals > 10e-16))/10000
            lower2 <- nonzero[1]
            upper2 <- nonzero[2]
            pval <- integrate(ifunc, lower2, upper2,...)$value
            } else {pval <- 0}
          }
          
        } else if (Delta<0){
          if (m==0){
            pval<- 1
          } else if (x==m){
            # x=m => BU ~ Beta(x+1,0)=1, so Pr[TL(2BU-1)>=Delta] = Pr[TL>=Delta]>=Pr[TL>=0]=1
            pval<-1
          } else {
            ifunc <- function(tt, delta = Delta) {
              pbeta(0.5 + delta / (2 * tt), x+1, m-x, lower.tail=FALSE) * dbeta(tt, m, n-m+1)
            }
            # do not integrate of the ifunc, unless there are at least 2 non-zero values (>1e-16)
            fvals <- ifunc(seq(0, 1, length.out = 10000))
            if(sum(fvals > 10e-16) > 1){
            nonzero <- range(which(fvals > 10e-16))/10000
            lower2 <- nonzero[1]
            upper2 <- nonzero[2]
            pval <- integrate(ifunc, lower2, upper2,...)$value
            } else {pval <- 0}
          }
          
        } 
      } else if (alt=="greater"){
        if (Delta==0){
          pval<- pbeta(0.5, x, m-x+1, lower.tail=TRUE)
        } else if (Delta>0){
          if (m==0){
            # m=0 => BL=0, TL=0 => Pr[TL(2BL-1) <=0] =1
            pval<-1
          } else if (x==0){
            # x=0, m>0 => BL=0 => Pr[TL(2BL-1)<=Delta]=Pr[-TL<=Delta]=1
            pval<-1
          } else {
            ifunc <- function(tt, delta = Delta) {
              pbeta(0.5 + delta / (2 * tt), x, m-x+1, lower.tail=TRUE) * dbeta(tt, m, n-m+1)
            }
            fvals <- ifunc(seq(0, 1, length.out = 10000))
            if(sum(fvals > 10e-16) > 1){
            nonzero <- range(which(fvals > 10e-16))/10000
            lower2 <- nonzero[1]
            upper2 <- nonzero[2]
            pval <- integrate(ifunc, lower2, upper2,...)$value
            } else {pval <- 0}
          }
          
        } else if (Delta<0){
          if (m==n){
            # m=n => TU~Beta(m+1,0)=1 => Pr[TU(2BL-1)<=Delta]=Pr[BL<=1/2 + Delta/2]
            pval<-pbeta(0.5 + Delta /2, x, m-x+1, lower.tail=TRUE)
          } else if (x==0){
            #x=0 => BL~Beta(0,m+1)=0 => Pr[TU(2BL-1)<=Delta]=Pr[-TU<=Delta]=Pr[TU>=-Delta]
            pval<- pbeta(-Delta, m+1,n-m, lower.tail=FALSE)
          } else {
            ifunc <- function(tt, delta = Delta) {
              pbeta(0.5 + delta / (2 * tt), x, m-x+1, lower.tail=TRUE) * dbeta(tt, m+1, n-m)
            }
            fvals <- ifunc(seq(0, 1, length.out = 10000))
            if(sum(fvals > 10e-16) > 1){
            nonzero <- range(which(fvals > 10e-16))/10000
            lower2 <- nonzero[1]
            upper2 <- nonzero[2]
            pval <- integrate(ifunc, lower2, upper2,...)$value
            } else {pval <- 0}
          }
          
        } 
        
      } else {
        stop("alt must be either 'less' or 'greater' ")
      }
      pval
    }
    
    
    
    rootFunc <- function(dd,XX,MM,NN, q, ALT) {
      pval.1sided(dd,X=XX, M=MM,N=NN, alt=ALT) - q
    }
    
    if (alternative=="two.sided" | alternative=="greater"){
      pL<- pval.1sided(nullparm,X=x,M=m,N=n, alt="greater")
      if (m==0){
        lower<- - qbeta(1-alpha,1,n)
      } else if (m==n & x==0){
        lower<- -1
      } else {
        lower<-  uniroot(rootFunc,interval = c(-1, 1),
                                          XX=x,MM=m,NN=n, 
                                          q = alpha, ALT="greater")$root
      } 
    }
    if (alternative=="two.sided" | alternative=="less"){
      pU<- pval.1sided(nullparm,X=x,M=m,N=n, alt="less")
      if (m==0){
        upper<- qbeta(1-alpha,1,n)
      } else if (m==n & x==m) {
        upper<- 1
      } else {
        upper<-  uniroot(rootFunc,interval = c(-1, 1),
                                          XX=x,MM=m,NN=n,
                                          q = alpha, ALT="less")$root
      }
      
    }
  } # end of nmc=0 algorithm
  
  ci <- c(lower = lower,
          upper = upper)
  attr(ci,which="conf.level")<-conf.level
  
  
  if (alternative=="two.sided"){
    p.value<- min(1,2*min(pL,pU))
  } else if (alternative=="greater"){
    p.value<- pL
  } else if (alternative=="less"){
    p.value<- pU
  }
  
  # create htest object
  dname <- paste0("n=",deparse(substitute(n)), " m=", deparse(substitute(m)),
                 " x=",deparse(substitute(x)))
  names(nullparm)<-"difference in proportions"
  METHOD<- "Exact McNemar Test (with central confidence intervals)"
  if (nmc>0) METHOD<-c(METHOD,paste0("Monte Carlo Implementation, nmc=",nmc,""))
  estimate<-c(x/n,(m-x)/n,x/n - (m-x)/n)
  names(estimate)<-c("x/n","(m-x)/n","difference")
  stat<-c(n,m,x)
  names(stat)<-c("n","m","x")
  #pvals<-c(pL0=pL0,pL=pL,pU0=pU0,pU=pU)
  #names(pvals)<-c("pL0","pL","pU0","pU")
  out<- list(statistic=stat,
             p.value=p.value, conf.int=ci,
             estimate=estimate,
             null.value=nullparm,
             alternative=alternative,
             method=METHOD,
             data.name=dname
  )
  class(out)<-"htest"
  out
}

#source("../mcnemarExactDP old.R")
#mcnemarExactDP_Edited(n=50,m=49,x=2)

#N<-50
#MMM<-49
#Y<-2
#mcnemarExactDP(Y,MMM,N)

#doChecks<-function(X,M,N){
#   # check: upper 95% limit should give p-value of 0.05 
#   out<-mcnemarExactDP(x=X,m=M,n=N, conf.level=.90)
#   upper<-mcnemarExactDP(x=X,m=M,n=N,alternative="less")$conf.int[2]
#   upCheck<-mcnemarExactDP(x=X,m=M,n=N,nullparm=upper,alternative="less")$p.value
#   # check: lower 95% limit should give p-value of 0.05 
#   lower<-mcnemarExactDP(x=X,m=M,n=N,alternative="greater")$conf.int[1]
#   downCheck<-mcnemarExactDP(x=X,m=M,n=N,nullparm=lower,alternative="greater")$p.value
#   list(out=out,upCheck=upCheck,downCheck=downCheck)
#}

# at edges, checks only work on one side.
#doChecks(0,20,20)
#doChecks(20,20,20)
#doChecks(0,0,20)