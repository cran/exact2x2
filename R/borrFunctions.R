



borrControl<-function(nAlphaGrid=10000,nThetaGrid=1000,  maxIter=0, digits=4, orderFunc=NULL){
  if (maxIter>100) warning("each iteration may take a long time, default maxIter=4")
  if (maxIter<0) stop("maxIter should be non-negative")
  if (!is.null(orderFunc) && !(orderFunc %in% c("AlphaGrid","ByRR"))) stop("orderFunc must be either NULL, 'AlphaGrid' or 'ByRR'  ")
  list(nAlphaGrid=nAlphaGrid,nThetaGrid=nThetaGrid,  maxIter=maxIter, digits=digits, orderFunc=orderFunc)
}


## Run on Biowulf (parallel processing computer) for all values n1=2,3,...,20 and n2=2,3,...,20
#setwd("H:/My Documents/methods/ZeroAdjusted2x2/R/borr Attempts")
#orderPreCalc<-readRDS("biowulf_matrix_list_complete.rds")
#NList<-seq(2,20)
#NCNT<-expand.grid(NList,NList)
#n1List<-NCNT[,1]
#n2List<-NCNT[,2]
#orderPreCalc<-list(n1List=NCNT[,1],n2List=NCNT[,2], tuningParm=0.025, orderList=orderPreCalc, 
#                   control=borrControl(nAlphaGrid = 10000, nThetaGrid=1000, maxIter=0, orderFunc="AlphaGrid"))
# check for best compression, it turns out that bzip2 is best here
#save(orderPreCalc, file="R/sysdata.rda")
#setwd("H:/My Documents/methods/PUBLISHED/2010 Fay Biostat (fisherci)/R/exact2x2")
#save(orderPreCalc, file="R/sysdata.rda",compress="bzip2")
#save(orderPreCalc, file="R/sysdata.rda",compress="xz")
#load("R/sysdata.rda")



# define functions
# getThreshold: gets threhold given tuningParm
# calcRejectProb: calculates maximum rejection prob given threshold
calcRejectProb <-
  function(p.ctrl,
           Threshold,
           p.trt = p.ctrl,
           n.trt,
           n.ctrl,
           max.uninf.ctrls = n.ctrl) {
    #Threshold should be the maximum number of ~INFECTED~ trt people if all controls infected
    p.reject <-
      (pbinom(Threshold[1], n.trt, p.trt)) * dbinom(n.ctrl, n.ctrl, p.ctrl)
    if (max.uninf.ctrls > 0) {
      for (i in 1:max.uninf.ctrls) {
        p.reject <-
          p.reject + (pbinom(Threshold[i + 1], n.trt, p.trt)) * dbinom(n.ctrl -
                                                                         i, n.ctrl, p.ctrl)
      }
    }
    return(p.reject)
  }
# getThreshold: gets threhold given tuningParm
getThreshold <-
  function(n.ctrl, n.trt,
           tuningParm = .025,
           nThetaGrid = 1000,
           max.uninf.ctrls = n.ctrl,
           forceConvex=TRUE) {
    
    
    
    
    p.seq <- seq(0, 1, length = nThetaGrid)
    thresholds <- rep(-1, max.uninf.ctrls + 1)
    
    #Choose rejection threshold for all controls infected
    
    for (i in 1:(max.uninf.ctrls + 1)) {
      continue <- TRUE
      
      while (continue == TRUE) {
        thresholds[i] <- thresholds[i] + 1
        error.next <-
          max(
            calcRejectProb(
              p.ctrl = p.seq,
              Threshold = thresholds[1:i],
              p.trt = p.seq,
              n.trt = n.trt,
              n.ctrl = n.ctrl,
              max.uninf.ctrls = i-1
            )
          )
        # ORIGINAL ALGORITHM did not check thresholds[i]=n.trt. Not really a problem 
        # untill you use this function to create the ordering for all values. Then it is a problem:
        #if (error.next > tuningParm | thresholds[i] >= n.trt | (i>1 && forceConvex && thresholds[i]>thresholds[i-1])) {
        #  
        # in the following && means we only check for greater than previous thresholds if i>1 and forceConvex
        if (error.next > tuningParm | thresholds[i] > n.trt | (i>1 && forceConvex && thresholds[i]>thresholds[i-1])) {
          continue <- FALSE
          thresholds[i] <- thresholds[i] - 1
        }
      }
    }
    # thresholds are given in eqn 1 of Gabriel et al for forceConxex=FALSE
    # threshold[1]=T_{alpha}^{N_C}, threshold[2]=T_{alpha}^{N_C-1}, etc.
    return(thresholds)
  }
# END OF Big FUNCTION DEFNITIONS



borrOrderingAlphaGrid <- function(n1,
                         n2,
                         tuningParm = .025,
                         controlborr=borrControl()) {
  
  nC<-n1
  nT<-n2

  
  nAlphaGrid<-controlborr$nAlphaGrid
  # we can figure out what the minAlpha is by solving for the p such that 
  # dbinom(0,nT,p)*dbinom(nC,nC,p)=(1-p)^nT * p^nC is maximized
  # taking the derivative, setting to zero, gives p=nC/(nC+nT)
  #   minAlpha<-controlborr$minAlpha
  minAlpha<- dbinom(0,nT,nC/(nC+nT))*dbinom(nC,nC,nC/(nC+nT))
  nThetaGrid<-controlborr$nThetaGrid
  maxIter<-controlborr$maxIter
  forceConvex<-TRUE
  
  threshold <- getThreshold(nC, nT, tuningParm, nThetaGrid=nThetaGrid, forceConvex=forceConvex)
  # binary.mat[i,j]=0 means reject if yC=nC-(i-1) and yT=j-1
  # here is the dimnames for the matrix, add it at the end
  # use 2 instead of 1, in case there are computer 'ties' at 1
  binary.mat <- matrix(2, nrow = nT + 1, ncol = nC + 1)
  for (i in 1:(nC + 1)) {
    if (threshold[i] >= 0) {
      binary.mat[0:(threshold[i] + 1), i] <- 0
    }
  }
  
  
  
  # take uniform sequence on log scale from 10^-7 to 1, and on the arithmetic scale from 0 to 1
  # then sort and take unique
  alpha.seq <-
    sort(unique(c(
      10 ^ seq(log10(minAlpha), 0, length = nAlphaGrid / 2),
      seq(minAlpha, 10 ^ 0, length = nAlphaGrid / 2)
    )))

  matNames<-list(paste0("x2=",0:nT),paste0("x1=",nC:0))
  
  getOrderMat<-function(alpha.seq){
    nAlpha <- length(alpha.seq)
    # get thresholds for different values of tuningParm
    thresholdMatrix <- matrix(NA, nrow = nAlpha, ncol = nC + 1)
    # we know that the threshold that goes with minAlpha is (0,-1,-1,...)
    thresholdMatrix[1,]<-c(0,rep(-1,nC))
    for (i in 2:nAlpha) {
      thresholdMatrix[i,] <- getThreshold(nC, nT, tuningParm = alpha.seq[i])
    }
    # estimate alpha values where each point enters rejection region 
    orderMat <- matrix(NA, nrow = nT + 1, ncol = nC + 1)
    for (j in 1:(nC + 1)) {
      for (i in 0:nT) {
        # Note old algorithm had: ==i not >=i 
        pick<- thresholdMatrix[, j] >= i
        if (any(pick)){
          orderMat[i + 1, j] <- min(alpha.seq[pick])
        } else {
          orderMat[i+1,j]<-1
          # previously had line below...but above is cleaner. If any(pick)=FALSE
          #               that means rejection should be at highest signf level
          #orderMat[i + 1, j] <- max(orderMat[i, j], orderMat[i + 1, j - 1])
        } 

      }
    }
    dimnames(orderMat)<-matNames
    return(orderMat)
  }
  
  orderMat<-getOrderMat(alpha.seq)
  # check for ties, if there are any, create a new alphaGrid with nAlphaGrid
  # equally spaced values between untied values
  if (forceConvex){
    # Do not count ties when x1=0 or x2=n2, since those will give tied p-values of 1
    ovec<-sort(as.vector(orderMat[-(nT+1),-(nC+1)]))
  } else {
    ovec<-sort(as.vector(orderMat))
    # will be ties at 1, ignore those
    # possible that may be extra ties at 1 that can be broken,
    # but practically that will not be a problem because those 
    # those values will have large and uninteresting p-values
    ovec<-ovec[ovec<1]
  }

  ties<-c(FALSE,diff(ovec)==0)
  iter<-0

  while (any(ties) & iter<maxIter){
    # get new alpha grid 
    utied<- unique(ovec[ties])
    ntied<-length(utied)
    newAlpha<-rep(NA,ntied*nAlphaGrid)
    cnt<-1
    for (i in 1:ntied){
      newAlpha[cnt:(cnt+nAlphaGrid-1)]<- seq(max(c(minAlpha,ovec[ovec<utied[i]])),
                                             min(c(1,ovec[ovec>utied[i]])), length=nAlphaGrid)
      cnt<-cnt+nAlphaGrid
    }
    newAlpha<-sort(unique(c(ovec,newAlpha,1)))
    # recalculate orderMat
    orderMat<-getOrderMat(newAlpha)
    # recheck for ties for next iteration
    ovec<-sort(as.vector(orderMat))  
    if (forceConvex){
      # Do not count ties when x1=0 or x2=n2, since those will give tied p-values of 1
      ovec<-sort(as.vector(orderMat[-(nT+1),-(nC+1)]))
    } else {
      # p-values might not be 1 if not convex. 
      ovec<-sort(as.vector(orderMat))
      ovec<-ovec[ovec<1]
    }
    ties<-c(FALSE,diff(ovec)==0)
    iter<-iter+1
  }
  if (iter==maxIter) warning("possible ties in the ordering function. Increasing nAlphaGrid or maxIter may break ties.")
  
  
  # values not in original rejection region that goes with the tuning parameter 
  # i.e., binary.mat=1
  # are ranked after all values in the original rejection region
  # see equation 5 of Gabriel et al 
  matNames<-list(paste0("x2=",0:nT),paste0("x1=",nC:0))
  mat <-
    matrix(
      rank(binary.mat + orderMat, ties.method = "first"),
      nrow = nT + 1,
      ncol = nC + 1, dimnames=matNames
    )
  # transpose and reverse rows so that the matrix matches borrPvals matrix
  mat<- t(mat)[(nC+1):1,]
  out<-list(rankMat=mat,alphaMat=t(orderMat)[(nC+1):1,])
  return(out)
  #return(mat)
}



borrOrderingByRR <-     function(n1,
                             n2,
                             tuningParm = .025,
                             controlborr=borrControl()) {
  
   ## This is Martha's function that orders by listing the rejection regions instead of 
   ## Searching through the alpha Grid
   ## It is much faster and more accurate (the AlphaGrid method is only as accurate as the 
   ## fineness of the grid) when n1 and n2 are small

  # switch to Martha's variable names
  nT<-n2
  nC<-n1
  alpha1<-tuningParm
  digits<- controlborr$digits
  nsteps<- controlborr$nThetaGrid
                            
  ## Define functions
  new.expand.grid <- function(input, reps) {
    expand.grid(replicate(reps, input, simplify = FALSE))
  }

  drop.concave <- function(x){
    p <- ncol(x)
    if(p==2) {filter <- x[,2]<=x[,1]
    } else {filter <- (x[,-1]<=x[,-p])%*%rep(1,p-1) ==(p-1)}
    x <- x[ filter, ]
    return(x) 
  }
  
  make.convex.grid <- function(nC,nT){
    
    count <- 2
    temp <- new.expand.grid(-1:nT,2)
    convex.set <- drop.concave(cbind(temp[,2],temp[,1]))
    
    while (count<=nC){
      count <- count+1
      section.size <- nrow(convex.set)
      convex.set <- do.call("rbind", replicate(nT+2,convex.set, simplify=FALSE))
      convex.set <- cbind(rep((-1):nT,each=section.size), convex.set )
      convex.set <- drop.concave(convex.set)
    }
    
    return(convex.set)
  }
  
  threshold <- getThreshold(n.ctrl=nC, n.trt=nT, tuningParm=alpha1)
  binary.mat <- matrix(1, nrow = nT + 1, ncol = nC + 1)
  
  for (i in 1:(nC + 1)) {
    if (threshold[i] >= 0) {
      binary.mat[0:(threshold[i] + 1), i] <- 0
    }
  }
  
  xmat <- make.convex.grid(nC=nC, nT=nT)
  
  p<-seq(0,1,length=nsteps)
  output <- apply(xmat, MARGIN=1, function(x) round(max(calcRejectProb(p,Threshold=x,n.trt=nT,n.ctrl=nC)),digits))
  
  output2 <- data.frame(cbind( xmat, output))
  
  #The following deletes rows, starting from the end of the matrix, in which the alpha level ("output") is not 
  #monotonically decreasing compared to the one that follows after
  squash <- T
  while (squash==TRUE){
    check1 <- c(output2[-1,]$output >output2[-nrow(output2),]$output,TRUE)
    if (min(check1)==0) {
      output2 <- output2[check1,]
    }else{squash <- FALSE}
  }
  
  
  order.mat <- matrix(NA, nrow = nT + 1, ncol = nC + 1)
  for (j in 1:(nC + 1)) {
    for (i in 0:nT) {
      order.mat[i + 1, j] <-   min(output2[output2[,j]>=i,]$output)
    }
  }
  
  mat <-    matrix( rank(binary.mat + order.mat, ties.method = "first"),
                    nrow = nT + 1,
                    ncol = nC + 1
  )
  
  matNames<-list(paste0("x2=",0:nT),paste0("x1=",nC:0))
  dimnames(mat)<-matNames
  dimnames(order.mat)<-matNames
  
  return(list(
    rankMat = t(mat)[(nC+1):1,],
    alphaMat = t(order.mat)[(nC+1):1,]
  ))
  
  
}




borrOrderingPreCalc<-function(n1,n2,tuningParm=0.025, orderPreCalc=orderPreCalc){
  n1List<-orderPreCalc$n1List
  n2List<-orderPreCalc$n2List
  if (orderPreCalc$tuningParm!=tuningParm) stop("tuningParm not equal to orderPreCalc$tuningParm")
  if (!(any(n1List==n1 & n2List==n2))) stop("must have n1 and n2 as  orderPreCalc$n1List[i] and orderPreCalc$n2List[i] for some i")
  m<- length(n1List)
  i<- (1:m)[n1List==n1 & n2List==n2]
  ordering<-orderPreCalc$orderList[[i]]
  
  if (!(nrow(ordering)==n1+1 & ncol(ordering)==n2+1)) stop("orderPreCalc does not match hard coded one. Need to edit borrOrderingPreCalc")


    dimnames(ordering)<-list(paste0("x1=",0:n1),paste0("x2=",0:n2))
  return(list(rankMat=ordering,alphaMat=NULL))
}

#borrOrderingPreCalc(2,6,tuningParm=0.025, orderPreCalc)


#borrOrderingPreCalc(3,6,orderPreCalc)
borrOrdering<-function (n1, n2, tuningParm = 0.025, controlborr = borrControl()) {
  orderFunc <- controlborr$orderFunc
  if (is.null(orderFunc)) {
    n1List <- orderPreCalc$n1List
    n2List <- orderPreCalc$n2List
    if (any(n1List == n1 & n2List == n2) & tuningParm == 
        orderPreCalc$tuningParm) {
      out <- borrOrderingPreCalc(n1, n2, tuningParm = tuningParm, 
                                 orderPreCalc)
    }
    else if (n1+n2<=16 & tuningParm != 
             orderPreCalc$tuningParm){
      out <- borrOrderingByRR(n1, n2, tuningParm = tuningParm, 
                              controlborr = controlborr)
    }
    else {
      out <- borrOrderingAlphaGrid(n1, n2, tuningParm = tuningParm, 
                                   controlborr = controlborr)
    }
  }
  else if (orderFunc == "AlphaGrid") {
    out <- borrOrderingAlphaGrid(n1, n2, tuningParm = tuningParm, 
                                 controlborr = controlborr)
  }
  else if (orderFunc == "ByRR") {
    out <- borrOrderingByRR(n1, n2, tuningParm = tuningParm, 
                            controlborr = controlborr)
  }
  return(out)
}


#b1<-borrOrdering(2,4,tuningParm=0.025, controlborr=
#               borrControl(nAlphaGrid = 10000, nThetaGrid=1000, maxIter=0, orderFunc="ByRR"))

#b2<-borrOrdering(4,4,tuningParm=0.025, controlborr=
#                   borrControl(nAlphaGrid = 10000, nThetaGrid=1000, maxIter=0, orderFunc="AlphaGrid"))
#b2$rankmat

#b3<-borrOrdering(2,4,tuningParm=0.025, controlborr=
#                   borrControl(nAlphaGrid = 10000, nThetaGrid=1000, maxIter=0, orderFunc=NULL))
#b3$rankMat
borrPreCalc <-
  function(NList=seq(2,20),
           tuningParm = 0.025,
           controlborr = borrControl()) {
    NCNT<-expand.grid(NList,NList)
    n1List<-NCNT[,1]
    n2List<-NCNT[,2]
    m <- length(n1List)
    orderList <- list()
    for (i in 1:m) {
        temp <-
          borrOrderingAlphaGrid(
            n1List[i],
            n2List[i],
            tuningParm = tuningParm,
            controlborr =controlborr
          )$rankMat
        dimnames(temp) <- NULL
        orderList <- c(orderList, list(temp))
        print(paste0("preCalc done: n1=", n1List[i], " n2=", n2List[i]))
    }
    

    return(
      list(
        n1List = n1List,
        n2List = n2List,
        tuningParm = tuningParm,
        orderList = orderList,
        control= controlborr
      )
    )
  }




borrTest<-function(x1,n1,x2,n2, tuningParm=0.025,    
                  parmtype = c("ratio", "difference","oddsratio"), nullparm = NULL, 
                  alternative = c("less", "greater","two.sided"),  
                  conf.int = TRUE, conf.level = 0.975, controlUC=ucControl(),
                  controlborr=borrControl(),...){
  # n1=nC and n2=nT

  if ( (n1>20 | n2>20 | tuningParm!=0.025) | (!is.null(controlborr$orderFunc) && controlborr$orderFunc=="AlphaGrid")) message("calculations may take a long time, see controlborr for options to increase speed, perhaps at the cost of accuracy")
  borrMat<-borrOrdering(n1,n2,tuningParm=tuningParm,
                        controlborr=controlborr)$rankMat
  Tstat.borrT <- function(X1, N1, X2, N2, delta0) {
    diag(borrMat[X1 + 1, X2 + 1])
  }
  alternative<-match.arg(alternative)
  parmtype<-match.arg(parmtype)
  if (alternative!="less") warning("alternative not less, borrTest only designed for alternative='less' ")
  
  # old versions of the program allowed forceConvex=FALSE
  # but it is no longer allowed
  controlborr$forceConvex<-TRUE
  if (controlborr$forceConvex & parmtype=="difference"){ ucMethod<-"user-fixed"
  } else { ucMethod<-"user"}
  out <-
    uncondExact2x2(x1,n1,x2,n2,
      method = ucMethod,
      parmtype = parmtype,
      alternative = alternative,
      Tfunc = Tstat.borrT,
      conf.int = conf.int,
      conf.level = conf.level, control=controlUC
    )
  if (parmtype=="difference"){
    description<-"Unconditional Exact Test on Difference in Proportions, method=borrT"
  } else if (parmtype=="ratio"){
    description<-"Unconditional Exact Test on Ratio of Proportions, method=borrT"
  } else if (parmtype=="oddsratio"){
    description<-"Unconditional Exact Test on Odds Ratio of Proportions, method=borrT"
  }
  out$method<-description
  out
}




borrPvals<-function(n1,n2, tuningParm=0.025,    
                  parmtype = c("ratio", "difference","oddsratio"), 
                  nullparm = NULL, alternative = c("less", "greater","two.sided"),  
                  conf.int = TRUE, conf.level = 0.975, 
                  controlUC=ucControl(), controlborr=borrControl(),...){
  # n1=nC and n2=nT
  # old versions of the program allowed forceConvex=FALSE
  # but it is no longer allowed
  controlborr$forceConvex<-TRUE
  borrMat<-borrOrdering(n1,n2,tuningParm=tuningParm,
                        controlborr=controlborr)$rankMat
  # 
  Tstat.borrT <- function(X1, N1, X2, N2, delta0) {
    diag(borrMat[X1 + 1,X2 + 1])
  }
  alternative<-match.arg(alternative)
  parmtype<-match.arg(parmtype)
  if (alternative!="less") warning("alternative not less, borrTest only designed for alternative='less' ")
  if (controlborr$forceConvex & parmtype=="difference"){ ucMethod<-"user-fixed"
  } else { ucMethod<-"user"}
  out <-
    uncondExact2x2Pvals(n1,n2,
                   method = ucMethod,
                   parmtype = parmtype,
                   alternative = alternative,
                   Tfunc = Tstat.borrT, control=controlUC)
  out
}

#library(exact2x2)
#borrTest(4,4,1,4)

#orderPreCalc<-borrPreCalc(n1List=2:20,
#                          n2List=2:20,
#                          tuningParm = 0.025,
#                          controlborr = borrControl(nAlphaGrid = 1000, nThetaGrid=10000, maxIter=0)) 
#save(orderPreCalc, file="R/sysdata.rda")


powerBorr<- function(n1,n2,p1,p2,alpha=0.025,...){
  pvals<-borrPvals(n1,n2,conf.int=FALSE,...)
  X1<- 0:n1
  X2<- 0:n2
  f1<- dbinom(X1,n1,p1)
  f2<- dbinom(X2,n2,p2)
  f<- matrix(rep(f1,n2+1)*rep(f2,each=n1+1),n1+1,n2+1)
  power<- sum(f[pvals<=alpha])
  power
}
