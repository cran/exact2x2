

power2gridRatio<-function(power2=3){
  denom<- 2^(power2-1)
  y<-0:denom
  c(y/denom,denom/rev(y)[-1])
}

power2gridDifference<-function(power2=3){
  denom<- 2^(power2-1)
  y<-0:denom
  c(y/denom-1,y[-1]/denom)
}

power2grid<-function(power2=3,from=10,to=1,dolog=TRUE){
  if (dolog){
    out<-exp(seq(from=log(from),to=log(to),length.out=1+2^power2))
  } else {
    out<-seq(from=from,to=to,length.out=1+2^power2)
  }
  out
}

unirootGrid<-function(func,power2=12, step.up=TRUE, pos.side=FALSE, print.steps=FALSE,power2grid=power2gridRatio,...){
  
  grid<- power2grid(power2)
  newf<-function(i){ 
    func(grid[i],...)      
  }
  uout<-uniroot.integer(newf,c(1,1+2^power2),step.power=power2-1,
                        print.steps = print.steps, step.up=step.up, pos.side=pos.side)
  if (uout$root==1){
    bound<-c(grid[1],grid[2])
  }  else if (uout$root==1+2^power2){
    bound<-c(grid[uout$root-1],grid[uout$root])
  } else {
    bound<-c(grid[uout$root-1],grid[uout$root+1])
  }
  list(iter=uout$iter,f.root= newf(uout$root), root=grid[uout$root], bound=bound)
}

#unirootGrid(function(x){ x-7.3}, print.steps=TRUE, power2grid=power2grid,pos.side = TRUE)



constrMLE.ratio<-function(X1,N1,X2,N2,rho0){
  # from Miettinen and Nurminen, 1985, see Stat Xact Procs 8, p. 298
  # I added in the limits as rho0=0 and rho0=Inf
  if (rho0==0){
    p1tilde<-(X1+X2)/(N1+X2)
    p2tilde<-rep(0,length(p1tilde))
  } else if (rho0==Inf){
    p2tilde<-(X1+X2)/(N2+X1)
    p1tilde<-rep(0,length(p2tilde))
  } else {
    A <- rho0 * (N1 + N2)
    B <- -1 * (rho0 * N2 + X2 + N1 + rho0 * X1)
    C <- X1 + X2
    p1tilde <- (-B - sqrt(B^2 - 4 * A * C))/(2 * A)
    p2tilde <- rho0 * p1tilde
  }
  # fix possible computer rounding errors
  # NEEDED, there is an error when: X1/N1=5/5 and X2/N2=0/7, rho0=1.000001e-06
  p1tilde<- pmin(1, p1tilde)
  p1tilde<- pmax(0, p1tilde)
  p2tilde<- pmin(1, p2tilde)
  p2tilde<- pmax(0, p2tilde)
  list(p1=p1tilde,p2=p2tilde)
}

#constrMLE.ratio(8,8,8,8,1)


constrMLE.difference<-function(X1,N1,X2,N2,delta0){
  ## Got this from Farrrington and Manning, Stat in Med 1990, 1447-1454
  ## they use Theta1-Theta2=delta instead of Theta2-Theta1=delta so switch
  ## for calculations, then switch back at the end
  N<-list(N1=N1,N2=N2)
  X<-list(X1=X1,X2=X2)
  X1 <- X$X2
  N1 <- N$N2
  X2 <- X$X1
  N2 <- N$N1
  P1hat <- X1/N1
  P2hat <- X2/N2
  # if (delta0==0){ p1d<-P1hat p2d<-P2hat } else { get MLEs (p1d and p2d)
  # from Appendix of Farrington and Manning, Stat in Med, 1990, 1447-1454
  theta <- N2/N1
  a <- 1 + theta
  s0 <- delta0
  b <- -1 * (1 + theta + P1hat + theta * P2hat + s0 * (theta + 
                                                         2))
  cc <- s0^2 + s0 * (2 * P1hat + theta + 1) + P1hat + theta * 
    P2hat
  d <- -P1hat * s0 * (1 + s0)
  v <- b^3/(3 * a)^3 - (b * cc)/(6 * a^2) + d/(2 * a)
  u <- sign(v) * (b^2/(3 * a)^2 - cc/(3 * a))^(1/2)
  temp <- v/u^3
  ## define 0/0=1 to avoid NaNs messing things up
  temp[v == 0 & u == 0] <- 1
  temp <- pmax(-1, temp)
  temp <- pmin(1, temp)
  w <- (1/3) * (pi + acos(temp))
  p1d <- 2 * u * cos(w) - b/(3 * a)
  p1d <- pmax(s0, p1d)
  p1d <- pmin(1, p1d)
  p2d <- p1d - s0
  p1d <- round(p1d, 8)
  p2d <- round(p2d, 8)
  # fix possible computer rounding errors (possibly not necessary)
  p1d<- pmin(1, p1d)
  p1d<- pmax(0, p1d)
  p2d<- pmin(1, p2d)
  p2d<- pmax(0, p2d)
  # switch back 
  list(p1=p2d,p2=p1d)
}


constrMLE.oddsratio<-function(X1,N1,X2,N2,delta0){
  # from Agresti and Min 2002, Biostatistics 3:379-386 use
  # P1=pihat_1(\theta_0)) delta0 = theta0 (in Agresti and Min) first
  # calculate the restricted MLE when OR=delta0 See Miettinen and
  # Nurminen, 1985, Stat in med 213-226 P1 = R0tilde (in MN) P2 = R1tilde
  # (in MN) N1=S0 N2=S1 X1+X2=c A = S0(OR-1)
  A <- N1 * (delta0 - 1)
  # B= S1*OR+S0 - c*(OR-1)
  B <- N2 * delta0 + N1 - (X1 + X2) * (delta0 - 1)
  C <- -(X1 + X2)
  
  ## to avoid computer errors (i.e., delta0=1+1e-15 giving very different values than delta0=1 )
  ## treat all as 1
  if (abs(delta0-1)<1e-10) {
    # when delta0=1 then A=0 and eqn is linear not quadratic no need for
    # quadratic formula
    delta0<-1
    P1 <- (X1 + X2)/(delta0 * N2 + N1)
  } else {
    P1 <- (-B + sqrt(B^2 - 4 * A * C))/(2 * A)
  }
  P2 <- (P1 * delta0)/(1 + P1 * (delta0 - 1))
  # fix possible computer rounding errors 
  P1<- pmin(1, P1)
  P1<- pmax(0, P1)
  P2<- pmin(1, P2)
  P2<- pmax(0, P2)
  list(p1=P1,p2=P2)
}




pickTstat<-function(method, parmtype="difference", tsmethod="central", 
                    alternative="two.sided"){
######## define the Tstat function
if (method=="FisherAdj"){
  # order function based on Fisher's exact mid-p at or=1
  Tstat<-function(X1,N1,X2,N2,delta0){
    phyper(X2,N2,N1,X2+X1) - 0.5*dhyper(X2,N2,N1,X1+X2)
  }   
} else  if (method == "score" & parmtype == "difference" & 
           tsmethod=="square" & alternative=="two.sided" ) {
  Tstat <- function(X1, N1, X2, N2, delta0) {
    temp<-constrMLE.difference(X1,N1,X2,N2,delta0)
    p1d<-temp$p1
    p2d<-temp$p2
    P1hat<-X1/N1
    P2hat<-X2/N2
    #  Now do Z round std.err and numerator to set fuzz to 0
    std.err <- round(((p1d * (1 - p1d))/N1 + (p2d * (1 - p2d))/N2)^(1/2), 
                     10)
    numerator <- round((P2hat - P1hat - delta0), 10)
    
    Z <- numerator/std.err
    ## set no effect to zero even as std.err goes to zero
    Z[numerator == 0 & std.err == 0] <- 0
    out <- Z^2
    out
  }
} else if (method == "score" & parmtype == "difference"){
  # Not tsmethod=="square & alternative=="two.sided"
  Tstat <- function(X1, N1, X2, N2, delta0) {
    temp<-constrMLE.difference(X1,N1,X2,N2,delta0)
    p1d<-temp$p1
    p2d<-temp$p2
    P1hat<-X1/N1
    P2hat<-X2/N2
    #  Now do Z round std.err and numerator to set fuzz to 0
    std.err <- round(((p1d * (1 - p1d))/N1 + (p2d * (1 - p2d))/N2)^(1/2), 
                     10)
    numerator <- round((P2hat - P1hat - delta0), 10)
    
    Z <- numerator/std.err
    ## set no effect to zero even as std.err goes to zero
    Z[numerator == 0 & std.err == 0] <- 0
    Z
  }
} else if (method == "simple" & parmtype =="difference" &
            tsmethod == "square" & alternative == "two.sided") {
  Tstat <- function(X1, N1, X2, N2, delta0) {
      stat <- (X2/N2) - (X1/N1) - delta0
      stat <- stat^2
      stat
  }
} else if (method == "simple" & parmtype =="difference"){
  # NOT tsmethod == "square" & alternative == "two.sided"
  Tstat <- function(X1, N1, X2, N2, delta0) {
    stat <- (X2/N2) - (X1/N1) - delta0
    stat
  }
} else if (method == "wald-pooled" & tsmethod == "square" & 
           alternative == "two.sided") {
  Tstat <- function(X1, N1, X2, N2, delta0) {
    P1hat <- X1/N1
    P2hat <- X2/N2
    Phat <- (X1 + X2)/(N1 + N2)
    numerator<-(P2hat - P1hat - delta0)
    denom<-(Phat * (1 - Phat) * (1/N1 + 1/N2))^(1/2)
    numerator<-round(numerator,10)
    denom<-round(denom,10)
    Z<- numerator/denom
    # if numerator and denom equal 0, set to 0
    Z[numerator==0 & denom==0]<-0
    Z <- Z^2
    Z
  }
} else if (method == "wald-pooled"){
  # NOT tsmethod == "square" & alternative == "two.sided"
  Tstat <- function(X1, N1, X2, N2, delta0) {
    P1hat <- X1/N1
    P2hat <- X2/N2
    Phat <- (X1 + X2)/(N1 + N2)
    numerator<-(P2hat - P1hat - delta0)
    denom<-(Phat * (1 - Phat) * (1/N1 + 1/N2))^(1/2)
    numerator<-round(numerator,10)
    denom<-round(denom,10)
    Z<- numerator/denom
    # if numerator and denom equal 0, set to 0
    Z[numerator==0 & denom==0]<-0
    Z
  }
} else if (method == "wald-unpooled" & tsmethod == "square" & 
           alternative == "two.sided") {
  Tstat <- function(X1, N1, X2, N2, delta0) {
    P1hat <- X1/N1
    P2hat <- X2/N2
    Phat <- (X1 + X2)/(N1 + N2)
    numerator<-(P2hat - P1hat - delta0)
    denom<- ((P1hat * (1 - P1hat))/N1+(P2hat * (1 - P2hat))/N2)^(1/2)
    numerator<-round(numerator,10)
    denom<-round(denom,10)
    Z<- numerator/denom
    # if numerator and denom equal 0, set to 0
    Z[numerator==0 & denom==0]<-0
    Z <- Z^2
    Z
  }
} else if (method == "wald-unpooled"){
  # NOT tsmethod == "square" & alternative == "two.sided"
  Tstat <- function(X1, N1, X2, N2, delta0) {
    P1hat <- X1/N1
    P2hat <- X2/N2
    Phat <- (X1 + X2)/(N1 + N2)
    numerator<-(P2hat - P1hat - delta0)
    denom<- ((P1hat * (1 - P1hat))/N1+(P2hat * (1 - P2hat))/N2)^(1/2)
    numerator<-round(numerator,10)
    denom<-round(denom,10)
    Z<- numerator/denom
    # if numerator and denom equal 0, set to 0
    Z[numerator==0 & denom==0]<-0
    Z
  }
} else if (method == "score" & parmtype == "oddsratio" &
           tsmethod == "square" & alternative == "two.sided") {
  Tstat <- function(X1, N1, X2, N2, delta0) {
    temp<-constrMLE.oddsratio(X1,N1,X2,N2,delta0)
    P1<-temp$p1
    P2<-temp$p2
    # use signed sqrt T in Agresti and Min, 2002 since we use
    # or=p2(1-p1)/p1(1-p2), switch first factor to n2,x2,p2, ect
    stat <- (N2 * (X2/N2 - P2))/(1/(N1 * P1 * (1 - P1)) + 1/(N2 * 
                                                               P2 * (1 - P2)))^(-0.5)
    # stat[(X1==0 & X2==0) | (X1==N1 & X2==N2)]<-0
    stat <- stat^2
    stat
  }
} else if (method == "score" & parmtype == "oddsratio"){
  # NOT tsmethod == "square" & alternative == "two.sided"
  Tstat <- function(X1, N1, X2, N2, delta0) {
    temp<-constrMLE.oddsratio(X1,N1,X2,N2,delta0)
    P1<-temp$p1
    P2<-temp$p2
    # use signed sqrt T in Agresti and Min, 2002 since we use
    # or=p2(1-p1)/p1(1-p2), switch first factor to n2,x2,p2, ect
    stat <- (N2 * (X2/N2 - P2))/(1/(N1 * P1 * (1 - P1)) + 1/(N2 * 
                                                               P2 * (1 - P2)))^(-0.5)
    # stat[(X1==0 & X2==0) | (X1==N1 & X2==N2)]<-0
    stat
  }
} else if (method == "score" & parmtype == "ratio" &
           tsmethod == "square" & alternative == "two.sided") {
  Tstat <- function(X1, N1, X2, N2, delta0) {
    # from Miettinen and Nurminen, 1985, see Stat Xact Procs 8, p. 298
    pi2hat <- X2/N2
    pi1hat <- X1/N1
    temp<-constrMLE.ratio(X1,N1,X2,N2,delta0)
    pi1tilde<-temp$p1
    pi2tilde<-temp$p2
    if (delta0==Inf){
      num<- -pi1hat
      denom<- rep(0,length(num))
    } else {
      num<-pi2hat - delta0 * pi1hat
      denom <- sqrt(
                 (pi2tilde * (1 - pi2tilde))/N2 + 
        (delta0^2 * pi1tilde * (1 - pi1tilde))/N1)
    }
    stat<- num/denom
    stat[num==0 & denom==0]<-0
    stat[X1==0 & X2==0]<-NA
    stat <- stat^2
    stat
  }
} else if (method == "score" & parmtype == "ratio"){
  # NOT tsmethod == "square" & alternative == "two.sided"
  Tstat <- function(X1, N1, X2, N2, delta0) {
    # from Miettinen and Nurminen, 1985, see Stat Xact Procs 8, p. 298
    pi2hat <- X2/N2
    pi1hat <- X1/N1
    temp<-constrMLE.ratio(X1,N1,X2,N2,delta0)
    pi1tilde<-temp$p1
    pi2tilde<-temp$p2
    if (delta0==Inf){
      num<- -pi1hat
      denom<- rep(0,length(num))
    } else {
      num<-pi2hat - delta0 * pi1hat
      denom <- sqrt(
        (pi2tilde * (1 - pi2tilde))/N2 + 
          (delta0^2 * pi1tilde * (1 - pi1tilde))/N1)
    }
    stat<- num/denom
    stat[num==0 & denom==0]<-0
    stat[X1==0 & X2==0]<-NA
    stat
  }
} else if (method == "simple" & parmtype =="ratio" &
           tsmethod == "square" & alternative == "two.sided") {
  Tstat <- function(X1, N1, X2, N2, delta0) {
    # originally, we had the following 2 lines, but 
    # we want the square to use log(theta2hat/(delta0*theta1hat))
    # so we can simplify
    # We put in the delta0 so when we square it, it makes sense
    # Similar to 'difference' equal theta2-theta1 - delta0
    # here we have: log(theta2) - log(theta1) - log(delta0)
    # 
    # Original (note ifelse not needed because x/0=Inf automatically
    #           and 0/0=NaN automatically):      
    #stat <- ifelse(X1 == 0, Inf, X2 * N1/(X1 * N2))
    #stat[X1 == 0 & X2 == 0] <- NA
    stat<- log(X2/N2) - log(X1/N1) - log(delta0) 
    stat <- stat^2
    stat
  }
} else if (method == "simple" & parmtype =="ratio"){
  # NOT tsmethod == "square" & alternative == "two.sided"
  Tstat <- function(X1, N1, X2, N2, delta0) {
    # originally, we had the following 2 lines, but 
    # we want the square to use log(theta2hat/(delta0*theta1hat))
    # so we can simplify
    # We put in the delta0 so when we square it, it makes sense
    # Similar to 'difference' equal theta2-theta1 - delta0
    # here we have: log(theta2) - log(theta1) - log(delta0)
    # 
    # Original (note ifelse not needed because x/0=Inf automatically
    #           and 0/0=NaN automatically):      
    #stat <- ifelse(X1 == 0, Inf, X2 * N1/(X1 * N2))
    #stat[X1 == 0 & X2 == 0] <- NA
    stat<- log(X2/N2) - log(X1/N1) - log(delta0) 
    stat
  }
} else if (method == "simple" & parmtype =="oddsratio" &
           tsmethod == "square" & alternative == "two.sided") {
  Tstat <- function(X1, N1, X2, N2, delta0) {
    # Note if T1=X1/N1 and T2=X2/N2 then: 
    #  T2(1-T1)/(T1(1-T2)) = X2(N1-X1)/(X1(N2-X2))
    stat<- log( X2*(N1-X1)/
                  (delta0*X1*(N2-X2))) 
    stat <- stat^2
    stat
  }
} else if (method == "simple" & parmtype =="oddsratio"){
  # NOT  tsmethod == "square" & alternative == "two.sided"
  Tstat <- function(X1, N1, X2, N2, delta0) {
    # Note if T1=X1/N1 and T2=X2/N2 then: 
    #  T2(1-T1)/(T1(1-T2)) = X2(N1-X1)/(X1(N2-X2))
    stat<- log( X2*(N1-X1)/
                  (delta0*X1*(N2-X2))) 
    stat
  }
}
}



calcTall<-function(Tstat, allx, n1, ally, n2, delta0=0, parmtype="difference", 
                   alternative="two.sided", tsmethod="central", EplusM=FALSE, tiebreak=FALSE){
  if (!EplusM & !tiebreak){ 
    return(Tstat(allx,n1,ally,n2,delta0))
  }  else if (tiebreak & tsmethod!="square"){
    ######### define function to break ties
    tie.f <- function(allx, n1, ally, n2, parmtype) {
      if (parmtype == "difference") {
        ### Break ties based on Z scores
        ### larger abs(Z scores) are treated as more extreme (further from 0) 
        theta1 <- allx/n1
        theta2 <- ally/n2
        std <- 1/(theta1 * (1 - theta1)/n1 + theta2 * (1 - theta2)/n2)^0.5
        new <- (theta2 - theta1)/(theta1 * (1 - theta1)/n1 + theta2 * 
                                    (1 - theta2)/n2)^0.5
        new[is.na(new)] <- 0  #when x1, x2 are 0 then std=0 and num=0 so set ratio to 0
      } else if (parmtype == "ratio") {
        # Ratio want large values to suggest high theta2/theta1
        # when x1=0 break ties by x2
        new <- ifelse(allx == 0 & ally > 0, ally, NA)
        # when x2=0 break ties by 1/x1 
        new <- ifelse(allx > 0 & ally == 0, 1/allx, new)
        # otherwise (unless x1=n1 and x2=n2) break ties by Z score on log(Ratio)
        new <- ifelse(allx > 0 & ally > 0, (log(ally) - log(n2) - log(allx) + 
                                              log(n1))/(1/allx - 1/n1 + 1/ally - 1/n2)^0.5, new)
        # x1=n1 and x2=n2 gives a Z score of 0/0=NaN, set to 0
        new[is.na(new)]<-0
      } else if (parmtype == "oddsratio") {
        # highest value of OR is suggested when x1=0 and x2=n2
        # but simple ties all x1=0 and x2=n2
        # break ties towards the point x=(0,n2)
        # lowest value of OR is suggested when x1=n1 and x2=0
        # break ties away from that point
        new <- ifelse(allx == 0 | allx==n1, ally/n2, NA)
        new <- ifelse(ally == n2 | ally==0,1-allx/n1, new)
        # Otherwise use Z score on log(OR)
        new <- ifelse((allx > 0 & allx < n1) & (ally > 0 & ally < n2), 
                      (log(ally) - log(n2 - ally) - log(allx) + log(n1 - allx))/(1/(allx) + 
                                                                                   1/(n1 - allx) + 1/(ally) + 1/(n2 - ally))^0.5, new)
        # set Z scores of 0/0 or Inf/Inf to 0
        new[is.na(new)]<-0 
      }
      new
    } 
    Talltemp<- Tstat(allx,n1,ally, n2, delta0)
    Talltiebreak<- tie.f(allx,n1,ally,n2,parmtype)
    # Use ranks to break ties, round to remove computer error
    R1<- rank(round(Talltemp,digits=10))
    R2<- rank(round(Talltiebreak,digits=10))
    # ranks are integers or values ending in 0.5
    # with min=1 and max=N =(n1+1)*(n2+1) (total sample size)
    # divide R2 by a factor so it is much less than 0.5
    # to break ties in R1
    out<- R1 + R2/(10*(n1+1)*(n2+1))
  } else if (tiebreak & tsmethod=="square"){
    stop("tiebreak=TRUE not supported with tsmethod='square' ")
  } 
    
  if (EplusM){
    if (tiebreak & tsmethod!="square"){
      Talltemp<- out
    } else {
      Talltemp<- Tstat(allx,n1,ally, n2, delta0)
    }
    # for x with no information, set to least extreme
    # set to -Inf when tsmethod=square
    if (alternative=="two.sided" & tsmethod=="square"){
      leastExtreme<- -Inf
    } else {
      leastExtreme<- Inf
    }
    if (parmtype=="ratio"){
        Talltemp[allx==0 & ally==0]<-  leastExtreme
    } else if (parmtype=="oddsratio"){
        Talltemp[(allx==0 & ally==0) | (allx==n1 & ally==n2)]<-  leastExtreme
    }
    out<-rep(NA,length(Talltemp))
    for (i in 1:length(Talltemp)){
      # get constrained MLE under null
      if (parmtype=="difference"){
        p1p2<-constrMLE.difference(allx[i],n1, ally[i], n2, delta0)
      } else if (parmtype=="ratio"){
        p1p2<-constrMLE.ratio(allx[i],n1, ally[i], n2, delta0)
      } else if (parmtype=="oddsratio"){
        p1p2<-constrMLE.oddsratio(allx[i],n1, ally[i], n2, delta0)
      } 
      
      if (alternative=="two.sided" & tsmethod=="square"){
        # for tsmethod='square' we want larger Talltemp[i] to give larger output
        # we want extreme to be larger T, so pval=sum(f[T>=t]) will give smaller
        # pval for larger T
        # so use 1-pval
        #
        I<-  Talltemp>=Talltemp[i]
        out[i]<- 1 - sum(dbinom(allx[I],n1,p1p2$p1)*dbinom(ally[I],n2,p1p2$p2) )
        
      } else {
        # get one-sided p-value so that order is the same direction
        # larger Talltemp[i] gives larger pval, and vise versa

        I<-  Talltemp<=Talltemp[i]
        out[i]<- sum(dbinom(allx[I],n1,p1p2$p1)*dbinom(ally[I],n2,p1p2$p2) )
      }
    }
  }
  out
}


getPrange<-function(x1,n1,x2,n2, parmtype, gamma){
  # get range for grid search over p1 and p2 when gamma>0
  if (gamma > 0) {
    bci1<-binom.test(x1,n1,conf.level=1-gamma/2)$conf.int
    bci2<-binom.test(x2,n2,conf.level=1-gamma/2)$conf.int
    p1min<-max(bci1[1],0)
    p1max<-min(bci1[2],1)
    p2min<-max(bci2[1],0)
    p2max<-min(bci2[2],1)            
  } else {
    p1min<-p2min<-0
    p1max<-p2max<-1
  }
  list(p1min=p1min, p1max=p1max, p2min=p2min, p2max=p2max)
}



getPI<-function(parmtype,DELTA, p1p2minmax, nPgrid){
    ### get range of PI1 values to search over range differs depending on the
    ### parmtype
    p1min<-p1p2minmax$p1min
    p2min<-p1p2minmax$p2min
    p1max<-p1p2minmax$p1max
    p2max<-p1p2minmax$p2max
    if (parmtype == "difference") {
        # p1 = p2-DELTA
        p1minNew<- max( p1min, p2min-DELTA)
        p1maxNew<- min( p1max, p2max-DELTA)
    } else if (parmtype == "ratio") {
        # p1 = p2/DELTA 
        p1minNew<- max( p1min, p2min/DELTA)
        p1maxNew<- min( p1max, p2max/DELTA)
    } else if (parmtype=="oddsratio"){
        # p1 = p2/(p2+DELTA-DELTA*p2)
        p1minNew<- max( p1min, p2min/(p2min+DELTA-DELTA*p2min))
        p1maxNew<- min( p1max, p2max/(p2max+DELTA-DELTA*p2max))
    }
    if (p1minNew==p1maxNew){
      PI1<- p1minNew
    } else if (p1minNew<p1maxNew){
      PI1 <- seq(p1minNew, p1maxNew, length = nPgrid)
    } else {
      PI1<-NA
    }
    
    if (parmtype == "difference") {
      # p1 = p2-DELTA
      PI2 <- PI1+DELTA
    } else if (parmtype == "ratio") {
      # p1 = p2/DELTA 
      PI2 <- PI1*DELTA
      # remove if PI1==0 and PI2==0
      keep<- !is.na(PI2/PI1) 
      PI1<-PI1[keep]
      PI2<-PI2[keep]
    } else if (parmtype=="oddsratio"){
      # p1 = p2/(p2+DELTA-DELTA*p2)
      PI2 <- DELTA*PI1/(1-PI1+DELTA*PI1)
      keep<- !is.na((PI2*(1-PI1))/(PI1*(1-PI2))) 
      PI1<-PI1[keep]
      PI2<-PI2[keep]
    }
    # if PI1=NA and PI2=NA, this means that Berger and Boos method excludes
    # all possible values
    list(PI1=PI1,PI2=PI2)
}
#pmm<-getPrange(3,10,5,12, parmtype="difference", gamma=10^-6, minPgridRatio=10^-6, 
#               maxPgridRatio=1-10^-6)

#PI<-getPI(parmtype="difference",DELTA=-1, pmm, nPgrid=100)
#PI
#PI
#PI$PI2-PI$PI1
#PI<-getPI(parmtype="ratio",DELTA=.004, p1min=0, p1max=1, p2min=0, p2max=1, nPgrid=5)
#PI
#PI$PI2/PI$PI1
#PI<-getPI(parmtype="oddsratio",DELTA=.4, p1min=0, p1max=1, p2min=0, p2max=1, nPgrid=5)
#PI
#(PI$PI2*(1-PI$PI1))/(PI$PI1*(1-PI$PI2))


getDeltaGrid<-function(parmtype, nCIgrid,maxPgridRatio=1-10^-6, minPgridRatio=10^-6){
  if (parmtype == "difference") {
    ds <- seq(-1, 1, length = nCIgrid)
  } else if (parmtype == "ratio") {
    delta.hi <- maxPgridRatio/minPgridRatio
    delta.low <- minPgridRatio/maxPgridRatio
    ds <- c(exp(seq(log(delta.low), log(delta.hi), length = nCIgrid)))
  } else if (parmtype == "oddsratio") {
    # make a grid to spread from 0 to Inf so that half are less than 1 and
    # half are greater
    delta.hi <- maxPgridRatio * (1 - minPgridRatio)/(minPgridRatio * (1 - maxPgridRatio))
    delta.low <- minPgridRatio * (1 - maxPgridRatio)/(maxPgridRatio * (1 - minPgridRatio))
    ds <- c(exp(seq(log(delta.low), log(delta.hi), length = nCIgrid)))
  }
  ds
}

#any(is.na(log(getDeltaGrid("ratio",1,100,1-10^-6,10^-6))))


######################### Define getPdQd function

getPDQD <- function(DELTA, Tstat, x1, n1, x2, n2, allx, ally, K1,K2, XX1, XX2,
                    parmtype, alternative, tsmethod, midp=FALSE, plotprobs=FALSE, EplusM=FALSE, tiebreak=FALSE, errbound=0,
                    p1p2minmax=NULL, nPgrid=100) {
  ####################################
  # first do Special Cases 
  ####################################
  # for DELTA=0 with ratio
  if (parmtype=="ratio" & DELTA==0){
    # DELTA=0 means pi2 must be zero, and pi1>0
    # if EplusM=FALSE and tiebreak=FALSE there may be lots of savings in 
    # not calculating Tall for x2>0, but for now do it the slow way 
    Tall<-calcTall(Tstat, allx, n1, ally, n2, delta0=DELTA, 
                   parmtype=parmtype, alternative=alternative, tsmethod=tsmethod, EplusM=EplusM, tiebreak=tiebreak)
    T0<-Tall[allx==x1 & ally==x2]
    # only pdf with positive values are when x2=0
    # also we condition on x!=(0,0)
    Tx2zero<- Tall[ally==0 & allx>0]
    if (all(Tx2zero<T0)){
      P<- 1
      Q<- 0
    } else if (all(Tx2zero>T0)){
      P<-0
      Q<-1
    } else {
      getPQdelta.eq.zero1<-function(p1){
        XX1<- allx[ally==0 & allx>0]
        # calculate the conditional pdf for Pr[X1=x | X1>0] 
        f1<- dbinom(XX1,n1,p1)/(1-dbinom(0,n1,p1))
        P <- sum(f1[Tx2zero < T0]) + sum((1 - (midp == TRUE) * 0.5) * 
                                       f1[Tx2zero == T0])
        Q <- sum(f1[Tx2zero > T0]) + sum((1 - (midp == TRUE) * 0.5) * 
                                       f1[Tx2zero == T0])
        list(P=P,Q=Q)
      }
      PI1 <- seq(p1p2minmax$p1min, p1p2minmax$p1max, length = nPgrid)
      PI1<- PI1[PI1>0]
      Pall<-Qall<-rep(NA,length(PI1))
      for (i in 1:length(PI1)){
        PQ<-getPQdelta.eq.zero1(PI1[i])
        Pall[i]<-PQ$P
        Qall[i]<-PQ$Q
      }
      P<-max(Pall)
      Q<-max(Qall)
    }
    return(list(P=P,Q=Q))
  } else if (parmtype=="ratio" & DELTA==Inf){
    # DELTA=Inf means pi2>0 and pi1=0
    # if EplusM=FALSE and tiebreak=FALSE there may be lots of savings in 
    # not calculating Tall for x2>0, but for now do it the slow way 
    Tall<-calcTall(Tstat, allx, n1, ally, n2, delta0=DELTA, 
                   parmtype=parmtype, alternative=alternative, tsmethod=tsmethod, EplusM=EplusM, tiebreak=tiebreak)
    T0<-Tall[allx==x1 & ally==x2]
    # only pdf with positive values are when x1=0
    # also we condition on x!=(0,0)
    Tx1zero<- Tall[allx==0 & ally>0]
    if (all(Tx1zero<T0)){
      P<- 1
      Q<- 0
    } else if (all(Tx1zero>T0)){
      P<-0
      Q<-1
    } else {
      getPQdelta.eq.zero2<-function(p2){
        XX2<- allx[allx==0 & ally>0]
        # calculate the conditional pdf for Pr[X1=x | X1>0] 
        f2<- dbinom(XX2,n2,p2)/(1-dbinom(0,n2,p2))
        P <- sum(f2[Tx1zero < T0]) + sum((1 - (midp == TRUE) * 0.5) * 
                                           f2[Tx1zero == T0])
        Q <- sum(f2[Tx1zero > T0]) + sum((1 - (midp == TRUE) * 0.5) * 
                                          f2[Tx1zero == T0])
        list(P=P,Q=Q)
      }
      PI2 <- seq(p1p2minmax$p2min, p1p2minmax$p2max, length = nPgrid)
      PI2<- PI2[PI2>0]
      Pall<-Qall<-rep(NA,length(PI2))
      for (i in 1:length(PI2)){
        PQ<-getPQdelta.eq.zero2(PI2[i])
        Pall[i]<-PQ$P
        Qall[i]<-PQ$Q
      }
      P<-max(Pall)
      Q<-max(Qall)
    }
    return(list(P=P,Q=Q))
  } else if (parmtype=="oddsratio" & (DELTA==0 | DELTA==Inf)){
    stop("Need to program getPDQD for parmtype='oddsratio' and DELTA=0 or Inf")
  }
  
  ###############################################################
  # End of Special Cases 
  ##############################################################
  if (errbound==0){
    Tall<-calcTall(Tstat, allx, n1, ally, n2, delta0=DELTA, 
                  parmtype=parmtype, alternative=alternative, tsmethod=tsmethod, EplusM=EplusM, tiebreak=tiebreak)
    T0<-Tall[allx==x1 & ally==x2]
  }
  ############## still inside get pdqd but defining getpq function
  getPQ <- function(pi1, pi2) {
    # see Stat Xact manual, e.g., StatXact8 Procs, p. 302 techinically all
    # 0<=pi1<=1 and 0<=pi2<=1, but in case there are computer numeric
    # issues, fix them now pi1[pi1>1]<-1 pi1[pi1<0]<-0 pi2[pi2>1]<-1
    # pi2[pi2<0]<-0 for speed do not use dbinom, use K1,K2, etc fx<-
    # rep(dbinom(0:n1,n1,pi1),n2+1) fy<-rep(dbinom(0:n2,n2,pi2),each=n1+1)
    if (errbound > 0) {
      ## dont search over very unlikely values
      XX1 <- max(0, qbinom(errbound/2, n1, pi1) - 1):min(n1, 
                                                         qbinom(1 - errbound/2, n1, pi1) + 1)
      XX2 <- max(0, qbinom(errbound/2, n2, pi2) - 1):min(n2, 
                                                         qbinom(1 - errbound/2, n2, pi2) + 1)
      # Aug 2, 2020: change to lchoose to avoid overflow
      K1 <- lchoose(n1, XX1)
      K2 <- lchoose(n2, XX2)
      allx <- rep(XX1, length(XX2))
      ally <- rep(XX2, each = length(XX1))
      Tall <- Tstat(allx, n1, ally, n2, DELTA)
      T0<- Tall[allx==x1 & ally==x2]
      
    }
    # Aug 2, 2020: fixed error 
    # here is an example that gave an incorrect p-value of 1
    # for versions <= 1.6.4.1
    # uncondExact2x2(x1=125, n1=1125,  x2=23, n2=30)
    # the true answer is a very small p-value
    #
    # to fix: changed K1 and K2 to lchoose to avoid overflow
    # so must change fx and fy also
    if (pi1==0 | pi1==1){
      fx.once<-  rep(0,length(XX1))
      if (pi1==0){
        fx.once[XX1==0]<- 1 
      } else if (pi1==1){
        fx.once[XX1==n1]<- 1 
      }
    } else {
      fx.once<- exp(K1 + XX1*log(pi1)+ (n1-XX1)*log(1-pi1))      
    }
    if (pi2==0 | pi2==1){
      fy.once<-  rep(0,length(XX2))
      if (pi2==0){
        fy.once[XX2==0]<- 1 
      } else if (pi2==1){
        fy.once[XX2==n2]<- 1 
      }
    } else {
      fy.once<- exp(K2 + XX2*log(pi2)+ (n2-XX2)*log(1-pi2))      
    }
    fx <- rep(fx.once, length(XX2))
    fy <- rep(fy.once, each = length(XX1))
    # end of Aug 2, 2020 fix
    # OLD CODE...Sometimes it would give K1<-choose() values of Inf so that f would have 
    # NaN values 
    #fx <- rep(K1 * pi1^XX1 * (1 - pi1)^(n1 - XX1), length(XX2))
    #fy <- rep(K2 * pi2^XX2 * (1 - pi2)^(n2 - XX2), each = length(XX1))
    
    f <- fx * fy
    ## sum(f) may be slightly less than 1 if errbound>0, standardize so sums to 1
    f<-f/sum(f)
    if (length(f)!=length(Tall)) stop("error in allx or ally")
    
    ## for ratio and odds ratio, there are some elements of the sample 
    ## space that give no information about the parameter
    ## e.g. x=(0,0) for ratio and additionally x=(n1,n2) for odds ratio
    ## We condition on data with information, and use the conditional 
    ## distribution
    if (parmtype=="ratio" | parmtype=="oddsratio"){
      if (parmtype=="ratio"){ 
        withInfo<-!(allx==0 & ally==0)
      } else {
        withInfo<- !((allx==0 & ally==0) | (allx==n1 & ally==n2))
      } 
      #if (any(is.na(Tall[withInfo]))) browser()
      allx<-allx[withInfo]
      ally<-ally[withInfo]
      Tall<-Tall[withInfo]
    
      # since x values without info have p-value=1, can never reject, and do not need to count them
      f<-f[withInfo]
    } 
    P<-Q<-NA
    if (any(is.na(Tall))){ 
      stop("NAs in Tstat function")}
    P <- sum(f[Tall < T0]) + sum((1 - (midp == TRUE) * 0.5) * 
                                     f[Tall == T0])
    Q <- sum(f[Tall > T0]) + sum((1 - (midp == TRUE) * 0.5) * 
                                     f[Tall == T0])
    Q <- ifelse(is.na(Q), 1, Q)
    P <- ifelse(is.na(P), 1, P)
    list(P = P, Q = Q)
  }
  ############# still inside get pdqd but end of getpq function
  
  ## PI1 = pi1: pi1 \in I(delta) see Stat Xact manual, e.g., StatXact8
  ## Procs, p. 302
  PI<-getPI(parmtype,DELTA, p1p2minmax, nPgrid)
  PI1<-PI$PI1
  PI2<-PI$PI2
  # if the gamma confidence intervals do not allow any of the DELTA values
  # then the PI1=PI2=NA
  # then set Pd and Qd to zero, p-values will add back gamma or gamma/2
  if (all(is.na(PI1)) | all(is.na(PI2))) return(list(Pd=0, Qd=0))
  
  Pall <- Qall <- rep(NA, length(PI1))
  for (i in 1:length(PI1)) {     
    PQ <- getPQ(PI1[i], PI2[i])  ###get probabilities over a grid of values of p1 and p2 defined by delta=p2-p1 (use dif as example)
    Pall[i] <- PQ$P  #### this is sum of prob less than observed
    Qall[i] <- PQ$Q  ##### sum of prob greater than observed
  }
  if (plotprobs){
    par(mfrow=c(2,1))
    plot(PI1, Qall,main=paste0("DELTA=",DELTA))
    plot(PI1, Pall,main=paste0("DELTA=",DELTA))
    par(mfrow=c(1,1))
  }
  Pd <- max(Pall)
  Qd <- max(Qall)
  out <- list(Pd = Pd, Qd = Qd)
  out
}
######################### End of getPDQD function


######## unconditional function method=='simpleTB' means simple tiebreak note
########## ngrid is the number of grid points for delta and nPgrid is the
########## number of probabilites to check over for each delta
ucControl<-function(nCIgrid = 500, errbound = 0, 
                    nPgrid = 100, power2=20, maxPgridRatio=1-10^-6, minPgridRatio=10^-6){
  
  # add checks if want
  if (minPgridRatio<=0) stop("minPgridRatio must be >0")
  if (maxPgridRatio>=1) stop("maxPgridRatio must be <1")
  if (nCIgrid<2) stop("nCIgrid must be >1")
  if (nPgrid<2) stop("nPgrid must be >1")
  list(nCIgrid=nCIgrid, errbound=errbound, nPgrid=nPgrid, power2=power2,
       maxPgridRatio=maxPgridRatio, minPgridRatio=minPgridRatio)
}



uncondExact2x2 <- function(x1, n1, x2, n2, 
    parmtype = c("difference", "ratio", "oddsratio"), nullparm = NULL, 
    alternative = c("two.sided","less", "greater"),  
    conf.int = FALSE, conf.level = 0.95, 
    method = c("FisherAdj", "simple", "score","wald-pooled", "wald-unpooled",  "user", "user-fixed"), 
    tsmethod = c("central","square"), midp = FALSE, 
    gamma = 0, EplusM=FALSE, tiebreak=FALSE,
    plotprobs = FALSE, control=ucControl(), Tfunc=NULL,...) {
  
    #dots<-match.call(expand.dots = TRUE)$"..."
    parmtype <- match.arg(parmtype)
    alternative <-match.arg(alternative)
    tsmethod <- match.arg(tsmethod)
    method <- match.arg(method)
    
    nCIgrid<-control$nCIgrid
    errbound<-control$errbound
    nPgrid<-control$nPgrid
    power2<-control$power2
    maxPgridRatio<-control$maxPgridRatio
    minPgridRatio<-control$minPgridRatio
    
    minparm <- switch(parmtype, difference = -1, ratio = 0, oddsratio = 0)
    maxparm <- switch(parmtype, difference = 1, ratio = Inf, oddsratio = Inf)

    if (is.null(nullparm)) {
      # set null hypothesis value of parameter if NULL to usual one
      nullparm <- switch(parmtype, difference = 0, ratio = 1, oddsratio = 1)
    }

    # when searching over range, cannot have Inf values, also Berger-Boos does not 
    # search over the whole range, so get p1 and p2 ranges
    p1p2minmax<-getPrange(x1,n1,x2,n2, parmtype, gamma)  
    
    

    ############################################
    # Check Conditions of Arguments
    ############################################
    if (plotprobs & conf.int) 
        stop("cannot have plotprobs=TRUE and conf.int=TRUE will produce too many plots")
    if (alternative!="two.sided" & tsmethod=="square") stop("when alternative='less' or 'greater', cannot use tsmethod='square'")
    if (tiebreak & tsmethod == "square") 
        stop("tie breaking function is designed to make sense only when tsmethod=central")
    if (gamma > 1 - conf.level) 
        stop("gamma is used for the Berger-Boos method and must be less than 1-conf.level, typically 10e-6")
    if (parmtype=="ratio" & (( (1/nullparm)< minPgridRatio) | (nullparm<minPgridRatio) ))
      stop("parmtype='ratio' and (nullparm<minPgridRatio) or (1/nullparm < minPgridRatio)")
    if (method=="user-fixed" & parmtype!="difference")
      stop("method='user-fixed' not allowed when parmtype!='difference' used method='user' ")
    if (method=="FisherAdj" & alternative=="two.sided" & tsmethod=="square")
      stop("tsmethod='square' does not work with method='FisherAdj', for tests like that use boschloo function")
    if (tiebreak & method!="simple") warning("tiebreak designed for method='simple', might not make sense for all methods")
    ######################### calculate variables that are needed inside functions
    ########################   but we only want to calculate once to save time
    ### To avoid clutter in the arguments of the functions, we do not put these variables as arguments
    

    allx <- rep(0:n1, n2 + 1)
    ally <- rep(0:n2, each = n1 + 1)

    
    
        
    ## for speed do the following only once, then use it instead of dbinom
    ## later
    if (errbound == 0) {
      XX1 <- 0:n1
      XX2 <- 0:n2
      # Aug 2, 2020: change to lchoose to avoid overflow
      # need to change the way fx and fy are calculated later 
      # to make it work out correctly
      K1 <- lchoose(n1, XX1)
      K2 <- lchoose(n2, XX2)
    }
    


    
    ### calculate estimate of parameter, and give names for output
    if (parmtype == "difference") {
      estimate <- (x2/n2) - (x1/n1)
      attr(estimate, "names") <- "p2-p1"
      attr(nullparm, "names") <- "p2-p1"
      description <- paste("Unconditional Exact Test on Difference in Proportions, method=", 
                           method)
    } else if (parmtype == "ratio") {
      if (x1 + x2 > 0) {
        estimate <- (x2/n2)/(x1/n1)
      } else {
        estimate <- NaN
      }
      attr(estimate, "names") <- "p2/p1"
      attr(nullparm, "names") <- "p2/p1"
      description <- paste("Unconditional Exact Test on Ratio of Proportions, method=", 
                           method)
    } else if (parmtype == "oddsratio") {
      p1 <- x1/n1
      p2 <- x2/n2
      if (x1 + x2 > 0) {
        estimate <- p2 * (1 - p1)/(p1 * (1 - p2))
      } else {
        estimate <- NaN
      }
      attr(estimate, "names") <- "p2(1-p1)/[p1(1-p2)]"
      attr(nullparm, "names") <- "p2(1-p1)/[p1(1-p2)]"
      description <- paste("Unconditional Exact Test on Odds Ratio, method=", 
                           method)
    }
    
    if (alternative=="two.sided" & tsmethod=="central"){
      description<-paste0(description,", central")
    } else if (alternative=="two.sided" & tsmethod=="square"){
      description<-paste0(description,", squared")
    } 
    if (EplusM){
      description<-paste0(description,", E+M")
    }
    if (tiebreak){
      description<-paste0(description,", with tie break")
      
    }
    ########################################### DEFINING FUNCTIONS
    
    
    
    ## Define Tstat Function
    if ( (method=="user" | method=="user-fixed") & is.null(Tfunc)){
      stop("when method='user' or 'user-fixed', you must supply Tfunc")
    } else if ( (method=="user" | method=="user-fixed") & !is.null(Tfunc)){
      ## We want Tstat<-Tfunc
      ## except it can be complicated if arguments are passed to Tfunc by ...
      ## first get list that is ... from call
      dots<-match.call(expand.dots = FALSE)$"..."
      if (length(dots)==0){
        Tstat<-Tfunc
      } else {
        ## Create a function to change the dots list into a character string as it was entered at ...         
        dotsAsChar<-function(dots){
          out<-","
          n<-length(dots)
          if (n==0){
            return("")
          } else if (n==1){
            out<-paste0(out, names(dots)[1],"=",dots[[1]]) 
          } else {
            for (i in 1:(n-1)){
              out<-paste0(out, names(dots)[i],"=",dots[[i]],",")
            }
            out<-paste0(out, names(dots)[n],"=",dots[[n]])
          }
          out
        }
        ## Create the cmd character string
        ## then evaluate it
        ## See p. 65 Venebles and Ripley 2000
        cmd<- paste("Tstat<-function(X1,N1,X2,N2,delta0){ 
                    Tfunc(X1,N1,X2,N2,delta0",dotsAsChar(dots),") }")
        eval(parse(text=cmd))
        # check that Tstat function is standard format
        if (alternative=="two.sided" & tsmethod=="square"){
          # want extremes to be larger than middle point
          Tmiddle<- Tstat(round(n1/2),n1,round(n2/2),n2)
          if (Tstat(0,n1,n2,n2)< Tmiddle | 
              Tstat(n1,n1,0,n2)< Tmiddle) warning("Tstat function appears non-standard. 'middle' (x1,x2) value gives T greater than at least one extreme (x1,x2) value")
        } else {
          if (Tstat(0,n1,n2,n2)<  Tstat(n1,n1,0,n2)) warning("Tstat function appears non-standard. Want T(x) at x=[x1,x2]=[0,n2] to be greater than T(x) at x=[n1,0]. Consider defining T(x) as -T(x).")
        }
    }
  } else {
      # non user supplied functions
      Tstat<-pickTstat(method=method, parmtype=parmtype, 
                       tsmethod=tsmethod, alternative=alternative)
  }
    


    
    ##  define getPdQd so we can set all the arguments except DELTA
    getPdQd <- function(DELTA){
         getPDQD(DELTA, Tstat=Tstat, x1=x1, n1=n1, x2=x2, n2=n2, 
                        allx=allx, ally=ally, K1=K1,K2=K2, XX1=XX1, XX2=XX2,
                        parmtype=parmtype, alternative=alternative, tsmethod=tsmethod, 
                        midp=midp, plotprobs=plotprobs, 
                        EplusM=EplusM, tiebreak=tiebreak, errbound=errbound,
                        p1p2minmax=p1p2minmax, nPgrid=nPgrid)
      
    }
  

                ######### define root function TO calculate CIs
        root.lower.f <- function(Delta, alpha) {
          if (parmtype!="difference" & Delta==0) return(-1)
          temp <- getPdQd(Delta)
          temp$Qd - alpha
        }
        root.upper.f <- function(Delta, alpha) {
          if (Delta==Inf) return(-1)
          temp <- getPdQd(Delta)
          temp$Pd - alpha
        }
        
        ######### define function TO calculate CI when the Barnard convexity conditions
        ######### hold and the ordering function (the Tstat function) does not change with delta0
        ######### This is a little faster than the other one
        conf.int.uniroot.f <- function(parmtype, conf.level, alternative) {
          alpha <- 1 - conf.level
          if (alternative == "two.sided") {
            alpha <- alpha/2
          }
          ci <- c(minparm, maxparm)
          if (alternative == "greater" | alternative == "two.sided") {
            # try uniroot, if it fails set $root value to minparm
            if (parmtype=="difference"){
              #temp <- tryCatch(uniroot(root.lower.f, c(minds, estimate), 
              #    alpha = alpha, maxiter = 500, extendInt="yes"), 
              #    error = function(e) list(root = minparm))
              if (sign(root.lower.f(-1, alpha=alpha))==sign(root.lower.f(1, alpha=alpha))) {
                temp<-list(root=-1)
              } else {
                temp<-unirootGrid(root.lower.f, power2, 
                                  power2grid=power2gridDifference, alpha=alpha)
              }
            } else {
              warning("unirootGrid may not work properly with parmtype!='difference' because of x values with no information, like x=(0,0)")
              temp<-unirootGrid(root.lower.f, power2, power2grid=power2gridRatio, alpha=alpha)
            }
            ci[1] <- temp$root
          }
          if (alternative == "less" | alternative == "two.sided") {
            # try uniroot, if it fails set $root value to maxparm
            if (parmtype=="difference"){
              #temp <- tryCatch(uniroot(root.upper.f, c(estimate, maxds), 
              #    alpha = alpha, maxiter = 500, extendInt="yes"), error = function(e) list(root = maxparm))
              if (sign(root.upper.f(-1, alpha=alpha))==sign(root.upper.f(1, alpha=alpha))) {
                temp<-list(root= 1)
              } else {
                temp <- unirootGrid(root.upper.f, power2=power2, power2grid = power2gridDifference, alpha = alpha)
              }
            } else {
              warning("unirootGrid may not work properly with parmtype!='difference' because of x values with no information, like x=(0,0)")
              temp <- unirootGrid(root.upper.f, power2=power2, power2grid= power2gridRatio, alpha = alpha)
            }
            ci[2] <- temp$root
          }
          ci
        }
        
        conf.int.nouniroot.f <- function(ds, conf.level, gamma, alternative, 
                                         tsmethod) {
          plower <- pupper <- rep(NA, length(ds))
          
          ## alpha is error, for one-sided we have 1-alpha=conf.level for
          ## two-sided we have 1-2*alpha=conf.level
          alpha <- 1 - conf.level - gamma
          if (alternative == "two.sided" & tsmethod == "central") 
            alpha <- alpha/2
          
          LOWER <- UPPER <- NA
          
          if (alternative == "greater" | (alternative == "two.sided" & tsmethod == 
                                          "central")) {
            for (i in 1:length(ds)) {
              out <- getPdQd(ds[i])
              plower[i] <- out$Qd
              if (out$Qd > alpha)
                break()
            }
            
            nomisslo <- !is.na(plower)
            i<- length(plower[nomisslo])
            if (i <= 1) {
              LOWER <- minparm
            }  else {
              # we know LOWER is between ds[i-1] and ds[i], where i=length(plower[nomisslo])
              # we could do a linear approximation...
              #LOWER <- switch(parmtype, 
              #  difference = approx(plower[nomisslo], ds[nomisslo], alpha)$y, 
              #  ratio = exp(approx(plower[nomisslo], log(ds[nomisslo]), alpha)$y), 
              #  oddsratio = exp(approx(plower[nomisslo], log(ds[nomisslo]), alpha)$y))
              if (parmtype=="difference"){
                p2g<-function(pow2){
                  power2grid(power2=pow2, from=ds[i-1], to=ds[i], dolog=FALSE)
                }
                temp<-unirootGrid(root.lower.f, power2, 
                                  power2grid=p2g, pos.side = TRUE, alpha=alpha)
              } else {
                p2g<-function(pow2){
                  power2grid(power2=pow2, from=ds[i-1], to=ds[i], dolog=TRUE)
                }
                temp<-unirootGrid(root.lower.f, power2, power2grid=p2g, pos.side = TRUE,  alpha=alpha)
              }
              LOWER<-temp$root
              
            }
          }
          
          if (alternative == "less" | (alternative == "two.sided" & tsmethod == 
                                       "central")) {
            # reverse order, start from largest and go down
            rds<-rev(ds)
            for (i in 1:length(rds)) {
              out <- getPdQd(rds[i])
              pupper[i] <- out$Pd
              if (out$Pd > alpha) break()
            }
            nomisshi <- !is.na(pupper)
            i<-length(pupper[nomisshi])
            if (i <= 1) {
              UPPER <- maxparm
              #warning(paste0("upper conf limit>", rds[1]," set to Inf. "))
            }  else {
              # we know UPPER is between rds[i-1] and rds[i]
              # we could do a linear approximation...
              #UPPER <- switch(parmtype, 
              #   difference = approx(pupper[nomisshi],ds[nomisshi], alpha)$y, 
              #   ratio = exp(approx(pupper[nomisshi], log(ds[nomisshi]), alpha)$y), 
              #   oddsratio = exp(approx(pupper[nomisshi],log(ds[nomisshi]), alpha)$y))
              if (parmtype=="difference"){
                p2g<-function(pow2){
                  power2grid(power2=pow2, from=rds[i], to=rds[i-1], dolog=FALSE)
                }
                temp<-unirootGrid(root.upper.f, power2, 
                                  power2grid=p2g, pos.side = TRUE, alpha=alpha)
              } else {
                p2g<-function(pow2){
                  power2grid(power2=pow2, from=rds[i], to=rds[i-1], dolog=TRUE)
                }
                temp<-unirootGrid(root.upper.f, power2, power2grid=p2g, pos.side = TRUE,  alpha=alpha)
              }
              UPPER<-temp$root
              
            }
            
            
            
  
            
          }
          if (alternative == "two.sided" & tsmethod == "square") {
            
            rds<-rev(ds)
            for (i in 1:length(rds)) {
              out <- getPdQd(rds[i])
              pupper[i] <- out$Qd
              
              if (pupper[i] > alpha) 
                  break()
            }
            nomisshi <- !is.na(pupper)
            i<- length(pupper[nomisshi])
            if (i<= 1) {
              UPPER <- maxparm
            } else {
              #UPPER <- switch(parmtype, 
              #  difference = approx(pupper[nomisshi], ds[nomisshi], alpha)$y, 
              #  ratio = exp(approx(pupper[nomisshi], log(ds[nomisshi]), alpha)$y), 
              #  oddsratio = exp(approx(pupper[nomisshi],log(ds[nomisshi]), alpha)$y))
              if (parmtype=="difference"){
                p2g<-function(pow2){
                  power2grid(power2=pow2, from=rds[i], to=rds[i-1], dolog=FALSE)
                }
                # use root.lower.f because want to use getPdQd()$Qd because it is 'square'

                temp<-unirootGrid(root.lower.f, power2, 
                                  power2grid=p2g, pos.side = TRUE, alpha=alpha)
              } else {
                p2g<-function(pow2){
                  power2grid(power2=pow2, from=rds[i], to=rds[i-1], dolog=TRUE)
                }
                temp<-unirootGrid(root.lower.f, power2, power2grid=p2g, pos.side = TRUE, alpha=alpha)
              }
              UPPER<-temp$root 
            }
            
            
            
            for (i in 1:length(ds)) {
              out <- getPdQd(ds[i])
              plower[i] <- out$Qd
              
              if (plower[i] > alpha) 
                 break()
            }
            nomisslo <- !is.na(plower)
            i<- length(plower[nomisslo])
            if (i<= 1) {
              LOWER <- minparm
            } else {
              #LOWER <- switch(parmtype, 
              #                difference = approx(plower[nomisslo],ds[nomisslo], alpha)$y, 
              #                ratio = exp(approx(plower[nomisslo], log(ds[nomisslo]), alpha)$y), 
              #                oddsratio = exp(approx(plower[nomisslo],log(ds[nomisslo]), alpha)$y))
              if (parmtype=="difference"){
                p2g<-function(pow2){
                  power2grid(power2=pow2, from=ds[i-1], to=ds[i], dolog=FALSE)
                }
                temp<-unirootGrid(root.lower.f, power2, 
                                  power2grid=p2g, pos.side = TRUE, alpha=alpha)
              } else {
                p2g<-function(pow2){
                  power2grid(power2=pow2, from=ds[i-1], to=ds[i], dolog=TRUE)
                }
                temp<-unirootGrid(root.lower.f, power2, power2grid=p2g, pos.side = TRUE,  alpha=alpha)
              }
              LOWER<-temp$root
            }
          }
          
          
          if (is.na(LOWER)) 
            LOWER <- minparm  
          if (is.na(UPPER)) 
            UPPER <- maxparm 
          
          ci <- c(LOWER, UPPER)
          ci
        }
        
        
        
        
        
        

        
############## End of Defining Functions ############################
        
        ################# calculating p values
        
        if ((parmtype=="ratio" & (x1+x2==0)) | 
            (parmtype=="oddsratio" & ((x1+x2==0) | (x1+x2==n1+n2)))){
            # Basically data have no information about parameter
            pval<-1
        } else if (alternative == "two.sided" & tsmethod == "central") {
            PQ0 <- getPdQd(nullparm)
            #pval <- min(1, 2 * min(PQ0$Pd + gamma/2, PQ0$Qd + gamma/2))
            pval <- min(1, 2 * min(PQ0$Pd + gamma, PQ0$Qd + gamma))
        } else if (alternative == "two.sided" & tsmethod == "square") {
            PQ0 <- getPdQd(nullparm)
            pval <- min(1,PQ0$Qd + gamma)
        } else if (alternative == "less") {
            PQ0 <- getPdQd(nullparm)
            pval <- min(1,PQ0$Pd + gamma)
        } else if (alternative == "greater") {
            PQ0 <- getPdQd(nullparm)
            pval <- min(1,PQ0$Qd + gamma)
        }
        
        
        ###############  find CI
        
        if (conf.int) {
            
            # with parmtype='ratio' or 'oddsratio' not sure we can show that the 
            # one-sided p-values are monotonic in the parameter
            if ((parmtype=="ratio" & (x1+x2==0)) | 
                (parmtype=="oddsratio" & ((x1+x2==0) | (x1+x2==n1+n2)))){
                  # Basically data have no information about parameter
                  ci <- c(minparm, maxparm)
            } else if ((parmtype=="difference") & (method == "simple" | 
              method=="user-fixed" | method=="FisherAdj") & (tsmethod == "central") & (gamma == 0)) {
                # fast algorith
                ci <- conf.int.uniroot.f(parmtype, conf.level, alternative)
            } else {
                # slower algorithm
                ds<-getDeltaGrid(parmtype, nCIgrid, maxPgridRatio, minPgridRatio)
                ci <- conf.int.nouniroot.f(ds, conf.level, gamma, alternative, 
                  tsmethod)
            }
            
        }
        
        
        ############### end of finding CI
        

    
    ############# preparing results
    if (conf.int == FALSE) 
        ci <- c(NA, NA)
    attr(ci, "conf.level") <- conf.level
    dname <- paste("x1/n1=(", x1, "/", n1, ") and x2/n2= (", x2, "/", n2, 
        ")", sep = "")

    statistic <- x1/n1
    names(statistic) <- "proportion 1"
    parameter <- x2/n2
    names(parameter) <- "proportion 2"
    out <- list(statistic = statistic, parameter = parameter, p.value = pval, 
        conf.int = ci, estimate = estimate, null.value = nullparm, alternative = alternative, 
        method = description, data.name = dname)
    class(out) <- "htest"
    out
}


uncondExact2x2Pvals<-function(n1,n2,...){
  allx<-rep(0:n1,n2+1)
  ally<-rep(0:n2, each=n1+1)
  pvals<-rep(NA, length(allx))
  for (i in 1:length(allx)){
      pvals[i]<-uncondExact2x2(allx[i],n1, ally[i], n2,...)$p.value    
  }
  out<- matrix(pvals,n1+1,n2+1)
  dimnames(out)<-  list(paste0("X1=",0:n1),paste0("X2=",0:n2))
  out
}

