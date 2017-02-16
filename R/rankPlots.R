createR<-function(Y){
    nr<-nrow(Y)
    nc<-ncol(Y)
    Y<-as.vector(Y)
    Y<-round(Y,7)
    rvec<-rank(Y,na.last="keep")
    rlen<-length(Y[!is.na(Y)])
    rstand<- (rvec-.5)/rlen
    R<-matrix(rstand,nr,nc)
    R
}


getColor<-function(r,COLORS=NULL,colorNA="black",ncolors=250){
    if (is.null(COLORS)){
        COLORS<-c(colorRampPalette(c("red","white"))(ncolors+3)[3:ncolors],"white",
               colorRampPalette(c("white","blue"))(ncolors+3)[2:(ncolors+1)])
     
    }  
    nC<-length(COLORS)
    coli<-round(nC*r)+1
    coli[coli<1]<-1
    coli[coli>nC]<-nC
    cout<-COLORS[coli]
    cout[is.na(cout)]<-colorNA
    cout
}

#getColor(0:10/10)
#plot(0:10/10,0:10/10,pch=16,col=getColor(0:10/10))


fillsquare<-function(i,j,r,n1,n2){
    I<-c(i-1,i-1,i,i,i-1)/(n1+1)
    J<-c(j-1,j,j,j-1,j-1)/(n2+1)
    if (is.na(r)){ COL<-"black"
    } else { COL<-getColor(r)
    }
    polygon(I,J,col=COL,border=NA)            
}

plotRankMat<-function(R,main=""){
    plot(c(0,1),c(0,1),type="n",main=main,xlab="", ylab="",axes=FALSE)
    n1<- nrow(R)-1
    n2<-ncol(R)-1
  
    at1<- 0:n1/(n1+1) + 0.5/(n1+1)
    X1<- 0:n1
    mtext(X1,side=1,at=at1,line=0,adj=0,cex=.8)
    at2<- 0:n2/(n2+1) + 0.5/(n2+1)
    X2<- 0:n2
    mtext(X2,side=2,at=at2,line=0,las=2,cex=.8)
    mtext(expression(X[1]),side=1,line=1.5,at=0.5)
    mtext(expression(X[2]),side=2,line=1.2,at=0.5,las=2)

    for (i in 1:(n1+1)){
        for (j in 1:(n2+1)){
            fillsquare(i,j,R[i,j],n1,n2)
        }
    }
}

plotT<-function(x,...){
  UseMethod("plotT")
}


plotT.function<-function(x,n1,n2,delta0=1, main="",...){
  Tstat<-x
  allx <- rep(0:n1, n2 + 1)
  ally <- rep(0:n2, each = n1 + 1)
  Tall<-Tstat(allx,n1,ally,n2,delta0,...)
  Tmat<-matrix(Tall,n1+1,n2+1)
  R<-createR(Tmat)
  plotRankMat(R,main=main)
}

plotT.numeric<-function(x,n1,n2,delta0=1, main="",...){
  Tall<-x
  if (length(Tall)!=(n1+1)*(n2+1)) stop("x vector must be of length (n1+1)*(n2+1)")
  Tmat<-matrix(Tall,n1+1,n2+1)
  R<-createR(Tmat)
  plotRankMat(R,main=main)
}

orderMat<-function(x,...){
  UseMethod("orderMat")
}

orderMat.numeric<-function(x,n1,n2,delta0,graphStyle=FALSE,...){
  Tall<-x
  if (length(Tall)!=(n1+1)*(n2+1)) stop("x vector must be of length (n1+1)*(n2+1)")
  out<- matrix(Tall,n1+1,n2+1)
  dimnames(out)<-  list(paste0("X1=",0:n1),paste0("X2=",0:n2))
  if (graphStyle) out<- out[(n1+1):1,]
  out
}

orderMat.function<-function(x,n1,n2,delta0,graphStyle=FALSE,...){
  Tfunc<-x
  X1<-rep(0:n1,n2+1)
  X2<-rep(0:n2,each=n1+1)
  out<- matrix(Tfunc(X1,n1,X2,n2,delta0,...),n1+1,n2+1)
  dimnames(out)<-  list(paste0("X1=",0:n1),paste0("X2=",0:n2))
  if (graphStyle) out<- out[(n1+1):1,]
  out
}

