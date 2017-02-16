binomMeld.test <-
  function(x1,n1,x2,n2,nullparm=NULL, 
           parmtype=c("difference","oddsratio","ratio"),
           conf.level=0.95, conf.int=TRUE,
           alternative=c("two.sided","less","greater"), midp=FALSE, nmc=0, eps=10^-8){
    
    ptype<-match.arg(parmtype)
    
    if (ptype=="difference"){
      g<-function(T1,T2){ T2-T1 }
      if (is.null(nullparm)) nullparm<-0
    } else if (ptype=="ratio"){
      g<-function(T1,T2){ T2/T1  }
      if (is.null(nullparm)) nullparm<-1
    } else if (ptype=="oddsratio"){
      g<-function(T1,T2){ T2*(1-T1)/( T1*(1-T2) ) }
      if (is.null(nullparm)) nullparm<-1
    }
    
    # use MLE, maybe other estimates may make sense but do not use them now
    estimate<-g(x1/n1,x2/n2)
    if (ptype=="difference"){
      names(estimate)<-"difference (p2-p1)"
    } else if (ptype=="ratio"){
      names(estimate)<-"ratio (p2/p1)"
    } else if (ptype=="oddsratio"){
      names(estimate)<-"odds ratio {p2(1-p1)}/{p1(1-p2)}"
    }  
    
    
    
    if (nmc==0){
      # calculate p.value and conf.int by numeric integration
      pvalCI<-binomMeldCalcInt(x1=x1,n1=n1,x2=x2,n2=n2,nullparm=nullparm, 
                               parmtype=parmtype,
                               conf.level=conf.level, conf.int=conf.int,
                               alternative=alternative, midp=midp, eps=eps)
    }  else {
      # calculate p.value and conf.int by Monte Carlo simulation, 
      # number of MC simulations=nmc
      pvalCI<-binomMeldCalcMC(x1=x1,n1=n1,x2=x2,n2=n2,nullparm=nullparm, 
                              parmtype=parmtype,
                              conf.level=conf.level, conf.int=conf.int,
                              alternative=alternative, midp=midp,nmc=nmc)
    }  
    
    dname<-paste("sample 1:(",x1,"/",n1,"), sample 2:(",x2,"/",n2,")",sep="")
    method<-paste("melded binomial test for",ptype)
    if (midp) method<-paste0(method,", mid-p version")
    if (nmc>0) method<-paste0(method,", by Monte Carlo, nmc=", nmc)
    stat<-x1/n1
    parm<-x2/n2
    names(stat) <- "proportion 1"
    names(parm) <- "proportion 2"
    names(nullparm)<-ptype
    alt<-match.arg(alternative)
    
    structure(list(statistic = stat, parameter = parm, 
                   p.value = pvalCI$p.value, 
                   conf.int = pvalCI$conf.int, estimate = estimate, null.value = nullparm, 
                   alternative = alt, method = method, 
                   data.name = dname), class = "htest")
    
  }




binomMeldCalcMC<-function(x1,n1,x2,n2,nullparm=NULL, 
                          parmtype=c("difference","oddsratio","ratio"),
                          conf.level=0.95, conf.int=TRUE,
                          alternative=c("two.sided","less","greater"),
                          midp=FALSE,nmc=10^6){
  ptype<-match.arg(parmtype)
  if (ptype=="difference"){
    g<-function(T1,T2){ T2-T1 }
    if (is.null(nullparm)) nullparm<-0
  } else if (ptype=="ratio"){
    g<-function(T1,T2){ T2/T1  }
    if (is.null(nullparm)) nullparm<-1
  } else if (ptype=="oddsratio"){
    g<-function(T1,T2){ T2*(1-T1)/( T1*(1-T2) ) }
    if (is.null(nullparm)) nullparm<-1
  }
  lowerLimit<- g(1,0)
  upperLimit<-g(0,1)
  
  alt<-match.arg(alternative)
  if (alt=="two.sided"){
    dolo<-dohi<-TRUE
    alpha<-(1-conf.level)/2  
  } else if (alt=="less"){
    ## alt=less so lower interval is lowest possible, do not calculate
    dolo<-FALSE
    lower<- lowerLimit
    dohi<-TRUE
    alpha<-1-conf.level
  } else if (alt=="greater"){
    # alt=greater so upper interval is highest possible, do not calculate
    dolo<-TRUE
    dohi<-FALSE
    upper<- upperLimit
    alpha<-1-conf.level
  } else stop("alternative must be 'two.sided', 'less', or 'greater' ")
  
  
  
  if (midp){
    U1<-rbinom(nmc,1,0.5)
    W1<- U1*rbeta(nmc,x1,n1-x1+1)+(1-U1)*rbeta(nmc,x1+1,n1-x1) 
    U2<-rbinom(nmc,1,0.5)
    W2<- U2*rbeta(nmc,x2,n2-x2+1)+(1-U2)*rbeta(nmc,x2+1,n2-x2)
    W1L<-W1U<-W1
    W2L<-W2U<-W2
  } else {
    # Not midp
    W1L<- rbeta(nmc,x1,n1-x1+1)
    W1U<- rbeta(nmc,x1+1,n1-x1) 
    W2L<- rbeta(nmc,x2,n2-x2+1)
    W2U<- rbeta(nmc,x2+1,n2-x2) 
  }
  
  if (alt=="greater" | alt=="two.sided"){
    gGreater<- g(W1U,W2L)
    # missing values should be defined as less than nullparm
    # since they should increase this p-value
    # e.g., all missing gives p=1, no information
    gGreater[is.na(gGreater)]<- lowerLimit
    pgr<- length( gGreater[gGreater<=nullparm])/nmc
  }
  if (alt=="less" | alt=="two.sided"){
    gLess<- g(W1L,W2U)
    # missing values should be defined as greater than nullparm
    # since they should increase this p-value
    # e.g., all missing gives p=1, no information
    gLess[is.na(gLess)]<- upperLimit
    pless<- length( gLess[gLess>=nullparm])/nmc
  }
  if (alt=="two.sided"){
    p.value<- min(1,2*pgr,2*pless)    
  } else if (alt=="less"){
    p.value<- pless
  } else if (alt=="greater"){
    p.value<- pgr
  }
  if (conf.int){
    lower<- lowerLimit
    upper<- upperLimit
    if (alt=="greater" | alt=="two.sided"){
      lower<- quantile(gGreater,probs=alpha)
    }
    if (alt=="less" | alt=="two.sided"){
      upper<- quantile(gLess, probs=1-alpha)
    }
  } else {
    lower<- NA
    upper<- NA
  }
  ci<-c(lower,upper)
  names(ci)<-NULL
  attr(ci,"conf.level")<-conf.level
  
  list(p.value=p.value, conf.int=ci)
}



binomMeldCalcInt <-
  function(x1,n1,x2,n2,nullparm=NULL, 
           parmtype=c("difference","oddsratio","ratio"),
           conf.level=0.95, conf.int=TRUE,
           alternative=c("two.sided","less","greater"), midp=FALSE, eps=10^-8){
    
    ptype<-match.arg(parmtype)
    if (ptype=="difference"){
      g<-function(T1,T2){ T2-T1 }
      if (is.null(nullparm)) nullparm<-0
    } else if (ptype=="ratio"){
      g<-function(T1,T2){ T2/T1  }
      if (is.null(nullparm)) nullparm<-1
    } else if (ptype=="oddsratio"){
      g<-function(T1,T2){ T2*(1-T1)/( T1*(1-T2) ) }
      if (is.null(nullparm)) nullparm<-1
    }
    
    lowerLimit<- g(1,0)
    upperLimit<-g(0,1)
    
    
    
    # Assume 
    # X1 ~ binom(n1,p1)
    # X2 ~ binom(n2,p2)
    #
    # with g(p1,p2)= p2-p1  (difference)
    #   or g(p1,p2)= p2/p1  (ratio)
    #   or g(p1,p2)= p2(1-p1)/(p1(1-p2))  (oddsratio)
    #
    #  Want to test with D=nullparm
    #    greater:   H0: g(p1,p2) <= D
    #               H1: g(p1,p2) > D             
    # or less:      H0: g(p1,p2) >= D
    #               H1: g(p1,p2) < D   
    #
    #
    # or two.sided, which has pvalue= min(1, 2*pg, 2*pl)
    #      where pg is p-value associated with "greater" alt Hyp
    #            pl is p-value associated with "less"  alt Hyp
    #   
    #    for greater (calculate lower CL) we use   
    # T1 ~ Beta(x1+1,n1-x1)
    # T2 ~ Beta(x2,n2-x2+1)
    #    and p-value is calculated under null: Pr[ g(T1,T2) <= D ]
    #
    #    for less (calculate upper CL) we use 
    # T1 ~ Beta(x1,n1-x1+1)
    # T2 ~ Beta(x2+1,n2-x2)
    #    and p-value is calculated under null: Pr[ g(T1,T2) >= D ]
    
    if (ptype=="difference"){ 
      # Pr[ g(T1,T2) >=D] = Pr[ T2-T1 >= D] = Pr[ T1 <= T2 - D]
      # = \int F1(t2 - D) f2(t2)  dt2
      # Use U1 and U2 for midp=TRUE
      #    T ~ Beta(x+U,n-x+1-U)
      # so U=0  means use lower confidence distribution
      #    U=1 means use upper CD
      # So if x2=0,U2=0 then f2 becomes a point mass at 0
      #     and if x2=n2, U2=1 then f2 becomes a point mass at 1
      # So we create a list that can do the appropriate integration
      # later: see Integrate function
      fLess<-list(
        f=function(t2,D,U1=0,U2=1){
          pbeta(t2-D,x1+U1,n1-x1+1-U1)*dbeta(t2,x2+U2,n2-x2+1-U2)
        },
        x20u20=function(D,x1,U1){ pbeta(-D,x1+U1,n1-x1+1-U1)},
        x2nu21=function(D,x1,U1){ pbeta(1-D,x1+U1,n1-x1+1-U1)},
        xnu1=function(D){ ifelse(0>=D,1,0)},
        x0u0=function(D){ ifelse(0>=D,1,0)}
      )
      
      # Pr[ g(T1,T2) <=D] = Pr[ T2-T1 <= D] = Pr[ T2 <= T1 + D]
      # = \int F2(t1 + D) f1(t1)  dt1
      fGreater<-list(
        f=function(t1,D,U1=1,U2=0){
          pbeta(t1+D,x2+U2,n2-x2+1-U2)*dbeta(t1,x1+U1,n1-x1+1-U1)
        },
        x10u10=function(D,x2,U2){ pbeta(D,x2+U2,n2-x2+1-U2)},
        x1nu11=function(D,x2,U2){ pbeta(1+D,x2+U2,n2-x2+1-U2)},
        xnu1=function(D){ ifelse(0<=D,1,0)},
        x0u0=function(D){ ifelse(0<=D,1,0)}
      )
    } else if (ptype=="ratio"){
      # Pr[ g(T1,T2) >=D] = Pr[ T2/T1 >= D] = Pr[ T1 <= T2/D]  
      # = \int F1(t2/D) f2(t2)  dt2
      fLess<-list(
        f=function(t2,D,U1=0,U2=1){
          pbeta(t2/D,x1+U1,n1-x1+1-U1)*dbeta(t2,x2+U2,n2-x2+1-U2)
        },
        x20u20=function(D,x1,U1){pbeta(0/D,x1+U1,n1-x1+1-U1)},
        x2nu21=function(D,x1,U1){ pbeta(1/D,x1+U1,n1-x1+1-U1)},
        xnu1=function(D){ ifelse(1>=D,1,0)},
        ## 0/0 gives no information, so always give p-value of 1
        x0u0=function(D){ 1 }
      )
      # Pr[ g(T1,T2) <=D] = Pr[ T2/T1 <= D] = Pr[ T2 <= T1*D]
      # = \int F2(t1 * D) f1(t1)  dt1
      fGreater<-list(
        f=function(t1,D, U1=1, U2=0){
          pbeta(t1*D,x2+U2,n2-x2+1-U2)*dbeta(t1,x1+U1,n1-x1+1-U1)
        },
        x10u10=function(D,x2,U2){pbeta(0,x2+U2,n2-x2+1-U2)},
        x1nu11=function(D,x2,U2){pbeta(D,x2+U2,n2-x2+1-U2)},
        xnu1=function(D){ ifelse(1<=D,1,0)},
        x0u0=function(D){ 1 }
      )  
    } else if (ptype=="oddsratio"){
      # Pr[ g(T1,T2) >=D] = Pr[ T2(1-T1)/T1(1-T2) >= D] 
      # = Pr[ T2(1-T1) >= T1(1-T2)D] = Pr[ T2 >= T1*(T2 + (1-T2)D) ]
      # = Pr[ T1 <= T2/{ T2 + (1-T2)D } ]  
      # = \int F1(t2/(t2+(1-t2)D) f2(t2)  dt2
      fLess<-list(
        f=function(t2,D, U1=0, U2=1){
          pbeta(t2/(t2+(1-t2)*D),x1+U1,n1-x1+1-U1)*dbeta(t2,x2+U2,n2-x2+1-U2)
        },
        x20u20=function(D,x1,U1){pbeta(0,x1+U1,n1-x1+1-U1)},
        x2nu21=function(D,x1,U1){ pbeta(1,x1+U1,n1-x1+1-U1)},
        ##  g(1,1)=0/0, no information
        xnu1=function(D){ 1 },
        x0u0=function(D){ 1 }
      )
      # Pr[ g(T1,T2) <=D] = Pr[ T2(1-T1)/T1(1-T2) <= D] 
      # = Pr[ T2(1-T1) <= T1(1-T2)D] = Pr[ T2(1-T1) + T1*T2*D <= T1*D  ]
      # = Pr[ T2( (1-T1) + T1*D) <= T1*D ]
      # = Pr[ T2 <= T1*D/{(1-T1) + T1*D} ]
      # = \int F2(t1*D/(1-t1+t1*D)) f1(t1)  dt1
      fGreater<-list(
        f=function(t1,D, U1=1, U2=0){
          pbeta(t1*D/(1-t1+t1*D),x2+U2,n2-x2+1-U2)*dbeta(t1,x1+U1,n1-x1+1-U1)
        },
        x10u10=function(D,x2,U2){ pbeta(0,x2+U2,n2-x2+1-U2)},
        x1nu11=function(D,x2,U2){pbeta(1,x2+U2,n2-x2+1-U2)},
        ##  g(1,1)=0/0, no information
        xnu1=function(D){ 1 },
        x0u0=function(D){ 1 }
      )
    }
    
    ## Create Integrate function, sometimes the integration is easy and is the 
    ## evaluation of a single function
    Integrate<-function(f, lo, hi, D, U1, U2){
      if (!is.null(f$xnu1) & x2==n2 & U2==1 & x1==n1 & U1==1){
        out<- f$xnu1(D)
      } else if (!is.null(f$x0u0) & x1==0 & U1==0 & x2==0 & U2==0){
        out<- f$x0u0(D)
      } else if (!is.null(f$x20u20) & x2==0 & U2==0){
        out<- f$x20u20(D, x1, U1)
      } else if (!is.null(f$x2nu21) & x2==n2 & U2==1){
        out<- f$x2nu21(D, x1, U1)
      } else if (!is.null(f$x10u10) & x1==0 & U1==0){
        out<- f$x10u10(D, x2, U2)
      } else if (!is.null(f$x1nu11) & x1==n1 & U1==1){
        out<- f$x1nu11(D, x2, U2)
      } else {
        out<-integrate(f$f, lo, hi, D=D, U1=U1, U2=U2)$value
      }
      list(value=out)
    }
    
    # p-value functions 
    pGreater<-function(delta){
      ## for the integrate function to work well, pick values that 
      ## make sense  in funcGreater
      ## recall funcGreater is 
      ##   pbeta2( W[t] )* dbeta1(t) 
      ##       where pbeta2(.)=pbeta(.,x2,n2-x2+1)
      ##             dbeta1(.)=dbeta(.,x1+1,n1-x1)
      ## and W[t] is different depending on the parmtype
      ##    difference: W[t] = t + D 
      ##       ratio:   W[t] = t*D
      ##   odds ratio:  W[t] = t*D/(1-t+t*D)
      ##
      getLimitsGr<-function(x1,n1,x2,n2,U1,U2,delta, eps){
        ## First, choose LowerInt=a and UpperInt=b so that \int_a^b  dbeta1(t) dt = 1- eps/2 
        LowerInt<-qbeta(eps/4,x1+U1,n1-x1+1-U1)
        UpperInt<-qbeta(1-eps/4,x1+U1,n1-x1+1-U1)
        ##  Second, choose a2 so that pbeta2( W[a2] )= eps/2
        ##       or   qbeta2( eps/2) = W[a2] 
        q<- qbeta(eps/2, x2+U2,n2-x2+1-U2)
        if (ptype=="difference"){
          a2<- q-delta
        } else if (ptype=="ratio"){
          a2<-q/delta
        } else if (ptype=="oddsratio"){
          # solve q = t*D/(1-t+t*D)   for t
          a2<- q/(delta+q-delta*q)
        }
        if (!is.na(a2)){
          LowerInt<-max(a2,LowerInt)
        } 
        list(lo=LowerInt,hi=UpperInt)
      }
      pout<-rep(0,length(delta))
      if (!midp){
        for (i in 1:length(delta)){
          limits<-getLimitsGr(x1,n1,x2,n2,U1=1,U2=0,delta=delta[i], eps=eps)
          if (limits$lo<limits$hi){
            pout[i]<-Integrate(fGreater,limits$lo,limits$hi,D=delta[i], U1=1, U2=0)$value
          }
        }
      } else {
        # midp=TRUE
        
        for (i in 1:length(delta)){
          limits00<-getLimitsGr(x1,n1,x2,n2,U1=0,U2=0,delta=delta[i], eps=eps)
          limits10<-getLimitsGr(x1,n1,x2,n2,U1=1,U2=0,delta=delta[i], eps=eps)
          limits01<-getLimitsGr(x1,n1,x2,n2,U1=0,U2=1,delta=delta[i], eps=eps)
          limits11<-getLimitsGr(x1,n1,x2,n2,U1=1,U2=1,delta=delta[i], eps=eps)
          
          pout[i]<- 0.25*Integrate(fGreater,limits00$lo,limits00$hi,
                                   D=delta[i], U1=0, U2=0)$value +
            0.25*Integrate(fGreater,limits10$lo,limits10$hi,
                           D=delta[i], U1=1, U2=0)$value +
            0.25*Integrate(fGreater,limits01$lo,limits01$hi,
                           D=delta[i], U1=0, U2=1)$value +
            0.25*Integrate(fGreater,limits11$lo,limits11$hi,
                           D=delta[i], U1=1, U2=1)$value 
          
        }
      }
      ## since we underestimate the integral (assuming perfect integration) by at most eps, 
      ## add back eps to get conservative p-value
      pout<-pout+eps
      pout
    }
    pLess<-function(delta){
      ## for the integrate function to work well, pick values that 
      ## make sense  in funcLess
      ## recall funcLess is 
      ##   pbeta1( W[t] )* dbeta2(t) 
      ##       where pbeta1(.)=pbeta(.,x1,n1-x1+1)
      ##             dbeta2(.)=dbeta(.,x2+1,n2-x2)
      ## and W[t] is different depending on the parmtype
      ##    difference: W[t] = t - D 
      ##       ratio:   W[t] = t/D
      ##   odds ratio:  W[t] = t/(t+(1-t)*D)
      ##
      getLimitsL<-function(x1,n1,x2,n2,U1,U2,delta, eps){
        ## First, choose LowerInt=a and UpperInt=b so that \int_a^b  dbeta2(t) dt = 1- eps/2 
        LowerInt<-qbeta(eps/4,x2+U2,n2-x2+1-U2)
        UpperInt<-qbeta(1-eps/4,x2+U2,n2-x2+1-U2)
        
        ##  Second, choose a1 so that pbeta1( W[a1] )= eps/2
        ##       or   qbeta1( eps/2) = W[a1] 
        q<- qbeta(eps/2, x1+U1,n1-x1+1-U1)
        if (ptype=="difference"){
          a2<- q+delta
        } else if (ptype=="ratio"){
          a2<-q*delta
        } else if (ptype=="oddsratio"){
          # solve q = t/(t+(1-t)*D)   for t
          a2<- q*delta/(1-q+delta*q)
        }
        if (!is.na(a2)) LowerInt<-max(a2,LowerInt)
        list(lo=LowerInt,hi=UpperInt)
      }
      pout<-rep(0,length(delta))
      if (!midp){
        for (i in 1:length(delta)){
          limits<-getLimitsL(x1,n1,x2,n2,U1=0,U2=1,delta=delta[i], eps=eps)
          if (limits$lo<limits$hi){
            pout[i]<-Integrate(fLess,limits$lo,limits$hi,D=delta[i], U1=0, U2=1)$value
          }
        }        
      } else {
        for (i in 1:length(delta)){
          
          limits00<-getLimitsL(x1,n1,x2,n2,U1=0,U2=0,delta=delta[i], eps=eps)
          limits10<-getLimitsL(x1,n1,x2,n2,U1=1,U2=0,delta=delta[i], eps=eps)
          limits01<-getLimitsL(x1,n1,x2,n2,U1=0,U2=1,delta=delta[i], eps=eps)
          limits11<-getLimitsL(x1,n1,x2,n2,U1=1,U2=1,delta=delta[i], eps=eps)
          pout[i]<- 0.25*Integrate(fLess,limits00$lo,limits00$hi,
                                   D=delta[i], U1=0, U2=0)$value +
            0.25*Integrate(fLess,limits10$lo,limits10$hi,
                           D=delta[i], U1=1, U2=0)$value +
            0.25*Integrate(fLess,limits01$lo,limits01$hi,
                           D=delta[i], U1=0, U2=1)$value +
            0.25*Integrate(fLess,limits11$lo,limits11$hi,
                           D=delta[i], U1=1, U2=1)$value 
        }
      }
      
      ## since we underestimate the integral (assuming perfect integration) by at most eps, 
      ## add back eps to get conservative p-value
      pout<-pout+eps
      pout
    }
    lower<-upper<-NA
    alt<-match.arg(alternative)
    if (alt=="two.sided"){
      dolo<-dohi<-TRUE
      alpha<-(1-conf.level)/2  
    } else if (alt=="less"){
      ## alt=less so lower interval is lowest possible, do not calculate
      dolo<-FALSE
      lower<- lowerLimit
      dohi<-TRUE
      alpha<-1-conf.level
    } else if (alt=="greater"){
      # alt=greater so upper interval is highest possible, do not calculate
      dolo<-TRUE
      dohi<-FALSE
      upper<- upperLimit
      alpha<-1-conf.level
    } else stop("alternative must be 'two.sided', 'less', or 'greater' ")
    
    
    
    
    
    
    if (dolo){
      ## Take care of special cases, when x2=0 T2 is point mass at 0
      ## when x1=n1 T1 is a point mass at 1
      if (!midp){
        if (x2==0 & x1<n1){
          if (conf.int) lower<- g( qbeta(1-alpha,x1+1,n1-x1), 0 )
          if (ptype=="difference"){ 
            pg<- 1- pbeta(-nullparm,x1+1,n1-x1)
          } else pg<-1
        } else if (x1==n1 & x2>0){
          if (conf.int) lower<- g( 1, qbeta(alpha,x2,n2-x2+1) )
          if (ptype=="difference"){
            ## Aug 22, 2014: Fixed following line
            ## WRONG: pg<- 1-pbeta(1+nullparm,x2,n2-x2+1)
            pg<- pbeta(1+nullparm,x2,n2-x2+1)
          } else if (ptype=="ratio"){
            pg<-pbeta(nullparm,x2,n2-x2+1)
          } else if (ptype=="oddsratio"){
            ## Aug 22, 2014: Fixed following line
            ## WRONG:pg<-pbeta(nullparm/(1-nullparm),x2,n2-x2+1)
            pg<-1
          }
        } else if (x2==0 & x1==n1){
          if (conf.int) lower<- g( 1, 0 )
          pg<-1
        } else {
          if (conf.int){ 
            rootfunc<-function(delta){
              pGreater(delta)-alpha
            }
            if (ptype!="difference"){  
              signRootLower<- sign(rootfunc(10^-100))
            } else {
              signRootLower<- sign(rootfunc(lowerLimit))
            }
            # signRootLower should be -1
            # because pGreater(delta) is increasing in delta
            if (signRootLower==1){
              # no root, so set lower confidence limit to lowerLimit
              lower<-lowerLimit
            } else {
              if (upperLimit==Inf){
                # uniroot cannot take Inf as an upper limit
                # find T such that rootfunc(T) has opposite sign as
                # rootfunc(lowerLimit)
                signRootUpper<- rootfunc(10^100)
                if (signRootLower==signRootUpper){
                  # no root, so set lower confidence limit to lowerLimit
                  lower<-lowerLimit
                } else {
                  for (i in 1:50){
                    upperLimit<-2^i
                    if (signRootLower!=sign(rootfunc(upperLimit))){
                      break()
                    }
                  }
                  
                  if (upperLimit==2^50){
                    warning("lower conf limit appears to be larger than 2^50=approx=10^16, set to 2^50")
                    lower<-2^50
                  } else {
                    lower<-uniroot(rootfunc,c(lowerLimit,upperLimit))$root
                  }
                }
              } else {
                # upperLimit<Inf
                lower<-uniroot(rootfunc,c(lowerLimit,upperLimit))$root
              }
            }
          }
          pg<-pGreater(nullparm)
          #reset upperLimit
          upperLimit<-g(0,1)
        }
      } else {
        # midp=TRUE
        if (conf.int){ 
          rootfunc<-function(delta){
            pGreater(delta)-alpha
          }
          if (ptype!="difference"){  
            signRootLower<- sign(rootfunc(10^-100))
          } else {
            signRootLower<- sign(rootfunc(lowerLimit))
          }
          
          if (signRootLower==1){
            lower<- lowerLimit
          } else {
            if (upperLimit==Inf){
              # uniroot cannot take Inf as an upper limit
              # find T such that rootfunc(T) has opposite sign as
              # rootfunc(lowerLimit)
              signRootUpper<- sign(rootfunc(10^100))
              if (signRootLower==signRootUpper){
                lower<-lowerLimit
              }  else {
                
                for (i in 1:50){
                  upperLimit<-2^i
                  if (signRootLower!=sign(rootfunc(upperLimit))){
                    break()
                  }
                }
              }
              
              if (upperLimit==2^50){
                warning("lower conf limit appears to be larger than 2^50=approx=10^15, set to 2^50")
                lower<-2^50
              } else {
                lower<-uniroot(rootfunc,c(lowerLimit,upperLimit))$root
              }
            } else {
              ## upperLimit<Inf
              lower<-uniroot(rootfunc,c(lowerLimit,upperLimit))$root
            }
          }
        }
        pg<-pGreater(nullparm)
        #reset upperLimit
        upperLimit<-g(0,1)
      }
    }     
    if (dohi){
      ## Take care of special cases, when x2=n2 T2 is point mass at 1
      ## when x1=0 T1 is a point mass at 0
      if (!midp){
        if (x1==0 & x2==n2){
          if (conf.int) upper<-g(0,1)
          pl<-1
        } else if (x1==0){
          if (conf.int) upper<- g(0,qbeta(1-alpha,x2+1,n2-x2))
          if (ptype=="difference"){
            pl<-1-pbeta(nullparm, x2+1,n2-x2)
          } else if (ptype=="ratio"){
            pl<- 1
          } else if (ptype=="oddsratio"){
            pl<-1
          }
        } else if (x2==n2){
          if (conf.int) upper<- g(qbeta(alpha,x1,n1-x1+1), 1)
          if (ptype=="difference"){
            ## Aug 22, 2014: Fixed following line
            ## WRONG:pl<-pbeta(1+nullparm, x1, n1-x1+1)
            pl<-pbeta(1-nullparm, x1, n1-x1+1)
          } else if (ptype=="ratio"){
            pl<-pbeta(1/nullparm,x1,n1-x1+1)
          } else if (ptype=="oddsratio"){
            pl<- 1
          }
        } else {
          if (conf.int){
            rootfunc<-function(delta){
              pLess(delta)-alpha
            }
            if (ptype!="difference"){  
              signRootLower<- sign(rootfunc(10^-100))
            } else {
              signRootLower<- sign(rootfunc(lowerLimit))
            }
            # signRootLower should be positive
            # because pLess(delta) is decreasing in delta
            if (signRootLower==-1){
              upper<-upperLimit
            } else {
              
              
              if (upperLimit==Inf){
                # uniroot cannot take Inf as an upper limit
                # find T such that rootfunc(T) has opposite sign as
                # rootfunc(lowerLimit)
                for (i in 1:50){
                  upperLimit<-2^i
                  if (signRootLower!=sign(rootfunc(upperLimit))){
                    break()
                  }
                }
              }
              if (upperLimit==2^50){
                warning("upper conf limit appears to be larger than 2^50=approx=10^15, set to Inf")
                upper<-Inf
              } else {
                upper<-uniroot(rootfunc,c(lowerLimit,upperLimit))$root
              }
            }
          }
          pl<-pLess(nullparm)
        }
      } else {
        # midp=TRUE
        if (conf.int){
          rootfunc<-function(delta){
            pLess(delta)-alpha
          }
          if (ptype!="difference"){  
            signRootLower<- sign(rootfunc(10^-100))
          } else {
            signRootLower<- sign(rootfunc(lowerLimit))
          }
          # signRootLower should be positive
          # because pLess(delta) is decreasing in delta
          if (signRootLower==-1){
            upper<-upperLimit
          } else {
            
            
            if (upperLimit==Inf){
              # uniroot cannot take Inf as an upper limit
              # find T such that rootfunc(T) has opposite sign as
              # rootfunc(lowerLimit)
              signRootUpper<- sign(rootfunc(10^100))
              if (signRootLower==signRootUpper){
                upper<-upperLimit
              } else {
                
                for (i in 1:50){
                  upperLimit<-2^i
                  if (signRootLower!=sign(rootfunc(upperLimit))){
                    break()
                  }
                }
                
                
                if (upperLimit==2^50){
                  warning("upper conf limit appears to be larger than 2^50=approx=10^15, set to Inf")
                  upper<-Inf
                } else {
                  upper<-uniroot(rootfunc,c(lowerLimit,upperLimit))$root
                }
              }
            } else {
              # upperLimit<Inf
              upper<-uniroot(rootfunc,c(lowerLimit,upperLimit))$root
              
            }
          }
        }
        pl<-pLess(nullparm)
      }
      
    }
    if (alt=="two.sided"){
      p.value<- min(1,2*pl,2*pg)
    } else if (alt=="less"){
      p.value<- pl
    } else if (alt=="greater"){
      p.value<- pg
    }
    if (conf.int){
      ci<-c(lower,upper)
    } else {
      ci<-c(NA,NA)
    }
    attr(ci,"conf.level")<-conf.level
    list(p.value=p.value,conf.int=ci)
  }

