boschloo<-function(x1, n1, x2, n2, alternative = c("two.sided", "less", "greater"), 
                   or = NULL, conf.int = FALSE, conf.level = 0.95, midp=FALSE, 
                   tsmethod=c("central","minlike"), control=ucControl()){
  
  alternative<-match.arg(alternative)
  tsmethod<-match.arg(tsmethod)
  
  
  if (midp){
    # when midp=TRUE, use the same ordering function for all
    # types of alternatives, but changes with delta0
    T<-function(X1,N1,X2,N2,delta0=1){
      m<-length(X1)
      pval<-rep(NA,m)
      for (i in 1:m){
        pval[i]<-exact2x2(matrix(c(X2[i],N2-X2[i],X1[i],N1-X1[i]),2,2),
                          or=delta0, conf.int=FALSE, 
                          midp=TRUE, alternative="less")$p.value
      }
      pval
    }
  } else if (alternative=="less"){
    # we want larger theta2hat to give higher T values
    T<-function(X1,N1,X2,N2,delta0=1){
      m<-length(X1)
      pval<-rep(NA,m)
      for (i in 1:m){
        pval[i]<-exact2x2(matrix(c(X2[i],N2-X2[i],X1[i],N1-X1[i]),2,2),
                          or=delta0, conf.int=FALSE, 
                          midp=FALSE, alternative="less")$p.value
      }
      pval
    }
  }  else if (alternative=="greater"){
    # we want larger theta2hat to give higher T values
    # but cannot use same T as for "less" because there are 
    # many ties at 1, so invert it
    T<-function(X1,N1,X2,N2,delta0=1){
      m<-length(X1)
      pval<-rep(NA,m)
      for (i in 1:m){
        pval[i]<-exact2x2(matrix(c(X2[i],N2-X2[i],X1[i],N1-X1[i]),2,2),
                          or=delta0, conf.int=FALSE, 
                          midp=FALSE, alternative="greater")$p.value
      }
      1-pval
    }
  }  else if (alternative=="two.sided" & tsmethod=="minlike"){
    T<-function(X1,N1,X2,N2,delta0=1){
      m<-length(X1)
      pval<-rep(NA,m)
      for (i in 1:m){
        pval[i]<-fisher.test(matrix(c(X2[i],N2-X2[i],X1[i],
                                      N1-X1[i]),2,2),or=delta0, conf.int=FALSE,
                             alternative="two.sided")$p.value
      }    
      1-pval
    }
    out<-uncondExact2x2(x1, n1, x2, n2, alternative = alternative, 
                        nullparm = or, parmtype = "oddsratio", conf.int = conf.int, 
                        conf.level = conf.level, midp=midp, method="user", 
                        tsmethod="square", control=control, Tfunc=T)
    out$method<-"Boschloo's test"
    return(out)
  } else if (alternative=="two.sided" & tsmethod=="central"){
    # recalcualte conf.level
    conf.level<- 1-(1-conf.level)/2
    # first do alternative="less"
    T<-function(X1,N1,X2,N2,delta0=1){
      m<-length(X1)
      pval<-rep(NA,m)
      for (i in 1:m){
        pval[i]<-exact2x2(matrix(c(X2[i],N2-X2[i],X1[i],N1-X1[i]),2,2),
                          or=delta0, conf.int=FALSE, 
                          midp=FALSE, alternative="less")$p.value
      }
      pval
    }
    outLess<-uncondExact2x2(x1, n1, x2, n2, alternative = "less", 
                            nullparm = or, parmtype = "oddsratio", conf.int = conf.int, 
                            conf.level = conf.level, midp=midp, method="user", control=control, Tfunc=T)
    # now do alternative="greater"
    T<-function(X1,N1,X2,N2,delta0=1){
      m<-length(X1)
      pval<-rep(NA,m)
      for (i in 1:m){
        pval[i]<-exact2x2(matrix(c(X2[i],N2-X2[i],X1[i],N1-X1[i]),2,2),
                          or=delta0, conf.int=FALSE, 
                          midp=FALSE, alternative="greater")$p.value
      }
      1-pval
    }
    outGreater<-uncondExact2x2(x1, n1, x2, n2, alternative = "greater", 
                               nullparm = or, parmtype = "oddsratio", conf.int = conf.int, 
                               conf.level = conf.level, midp=midp, method="user", control=control, Tfunc=T)
    out<-outLess
    out$p.value<- min(c(1,2*outLess$p.value,2*outGreater$p.value))
    out$conf.int<-c(outGreater$conf.int[1],outLess$conf.int[2])
    out$method<-"Boschloo's test"
    return(out)
  }
  uncondMethod<-"user"
  out<-uncondExact2x2(x1, n1, x2, n2, alternative = alternative, 
                      nullparm = or, parmtype = "oddsratio", conf.int = conf.int, 
                      conf.level = conf.level, midp=midp, method=uncondMethod, control=control, Tfunc=T)
  out$method<-"Boschloo's test"
  if (tsmethod=="central" & alternative=="two.sided") out$method<-"Central Boschloo's test"
  out
}

