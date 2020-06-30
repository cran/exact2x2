powerPaired2x2 <-
  function(pb,pc,npairs,sig.level=0.05,alternative=c("two.sided","one.sided"),strict=FALSE,
           nullOddsRatio=1,errbound=10^-6,...){
    prob.reject<-0
    dots<-match.call(expand.dots=FALSE)$"..."
    if (any(names(dots)=="or")) stop("use 'nullOddsRatio' instead of 'or' ")
    dotsDescription<- paste(names(dots),dots,sep="=", collapse = ",")
    
    if (pb>1 | pc>1 | pb<0 | pc<0) stop("pb and pc should be between 0 and 1")
    #TSMETHOD<-tsmethod
    eps<-errbound
    alt<-match.arg(alternative)
    ## SIG.LEVEL and ALT are values that will be used in calculation
    ## not output
    SIG.LEVEL<-sig.level
    ALT<-alt
    if (ALT=="less" | ALT=="greater"){
      stop("alternative must be either 'one.sided' or 'two-sided', the one-sided direction is defined by pb and pc to maximize power")
    }
    if (!strict & ALT=="two.sided"){
      # for two.sided (strict==FALSE) you want to calculate power on one-sided p-values at alpha/2
      # so you do not count rejections in the wrong direction
      SIG.LEVEL<- SIG.LEVEL/2
      ALT<-"one.sided"
    }
    if (ALT=="one.sided"){
      if (pb<pc){ ALT<-"less"
      } else { ALT<-"greater" } 
    }
    
    
    
    
    if (errbound==0){
      blow<-0
      bhigh<-npairs
      clow<-1
      chigh<- npairs
    } else {
      blow<-max(0,qbinom(eps/2,npairs,pb)-1)
      bhigh<-min(npairs,qbinom(1-eps/2,npairs,pb)+1)
      clow<-max(0,qbinom(eps/2,npairs,pc)-1)
      chigh<-min(npairs,qbinom(1-eps/2,npairs,pc)+1)
    }
    
    
    for (B in blow:bhigh){
      for (C in clow:chigh){
        # it would be slightly faster not put data in matrix form, but to keep from making programming errors,
        # just call exact2x2
        # it does not matter how the values that are not B or C are distributed
        # so just put them all in the A cell and let the D cell be 0
        notBC<- npairs- B-C 
        if (notBC<0) next()
        x<-matrix(c(notBC,B,C,0),byrow=TRUE,2,2)
        pval<-exact2x2(x,or=nullOddsRatio, alternative=ALT, conf.int=FALSE, paired=TRUE, plot=FALSE,...)$p.value
        #pval<-0
        if (pval<=SIG.LEVEL){ 
          prob.reject<-prob.reject+ dmultinom(c(notBC,B,C),size=npairs, prob=c(1-pb-pc,pb,pc))
        }
      }
    }
    
    ## create pretty output using power.htest class
    if (!strict){
      METHOD<-"Power for McNemar's Exact Test"
    } else if (strict){
      METHOD<-"Power for McNemar's Exact Test, including power to reject in the wrong direction"
    }
    
    if (!is.null(dots)) METHOD<-paste0(METHOD,", using (",dotsDescription,")")
    
    #if (approx) METHOD<- paste("Approximate",METHOD)
    
    ## NOTE: output original sig.level NOT SIG.LEVEL
    ##    e.g., two.sided sig.level=0.05 should be calculated (when strict=FALSE)
    ##    at level SIG.LEVEL=0.025 
    
    NOTE<-paste("errbound=",errbound)
    if (pb+pc==1){
      NOTE<-paste(NOTE," no non-informative pairs: pa+pd=0")
    }
    
    output<-list(power=prob.reject, npairs=npairs, 
                 pb=pb,pc=pc, sig.level = sig.level, 
                 alternative =alt, nullOddsRatio=nullOddsRatio, note =NOTE, 
                 method = METHOD)
    
    structure(output,class = "power.htest")
  }

#library(exact2x2)
#t0<-proc.time()
#powerPaired2x2(.5,.3,npairs=20,tsmethod=NULL)

#,tsmethod="central", midp=FALSE)
#t1<-proc.time()
#powerPaired2x2(.4,.3,npairs=500,errbound=0)
#t2<-proc.time()