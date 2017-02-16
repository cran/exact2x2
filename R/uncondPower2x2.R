uncondPower2x2<-function(n1,n2, theta1, theta2, alpha, ...){
  pval<-as.vector( uncondExact2x2Pvals(n1,n2,...) )
  f1<- dbinom(0:n1,n1, theta1)
  f2<- dbinom(0:n2, n2, theta2)
  f<- rep(f1,n2+1)*rep(f2, each=n1+1)
  
  power<- sum(f[pval<=alpha])
  power
}





Power2x2<-function(n1,n2, theta1, theta2, alpha, pvalFunc, ...){
  allx<-rep(0:n1,n2+1)
  ally<-rep(0:n2, each=n1+1)
  pvals<-rep(NA, length(allx))
  for (i in 1:length(allx)){
    pvals[i]<-pvalFunc(allx[i],n1, ally[i], n2)    
  }
  f1<- dbinom(0:n1,n1, theta1)
  f2<- dbinom(0:n2, n2, theta2)
  f<- rep(f1,n2+1)*rep(f2, each=n1+1)
  
  power<- sum(f[pvals<=alpha])
  power
}

SS2x2<-function(theta1, theta2, alpha, pvalFunc, power=0.90, n1start=10, increaseby=1, n2.over.n1=1,  
                maxiter=50, printSteps=TRUE, ...){
  n1<-n1start
  n2<- round(n2.over.n1*n1)
  pow<- Power2x2(n1,n2, theta1,theta2, alpha, pvalFunc, ...)
  if (printSteps){
    message("n1=",n1,"n2=",n2,"power=",pow)
  }
  if (pow>power){
    warning("power calculated at n1start is greater than 'power', there may be a smaller sample size")
  } else {
    for (i in 1:maxiter){
      n1<- n1+increaseby
      n2<- round(n2.over.n1*n1)
      pow<- Power2x2(n1,n2, theta1,theta2, alpha, pvalFunc, ...)
      if (printSteps){
        message("n1=",n1,"n2=",n2,"power=",pow)
      }
      
      if (pow>power){
        break()
      } 
      if (i==maxiter) warning("reached maxiter without achieving desired power")
    }
    
  }
  out<-list(power=pow,theta1=theta1,theta2=theta2,
            n1=n1,n2=n2,
            sig.level=alpha, method="rejection determined by pvalFunc")
  class(out)<-"power.htest"
  out
}

