powerUncondExact2x2<-function(n1,n2,theta1,theta2,
    alpha=0.05,...){
    ## Mike wrote July 13, 2016
    ## to do: perhaps set alpha=NULL and go to .05 if two.sided
    ## and .025 if one-sided

    ## perhaps add errbound for increasing speed
    ##  for large sample sizes 
    X1<- 0:n1
    X2<- 0:n2
    f1<- dbinom(X1,n1,theta1)
    f2<- dbinom(X2,n2,theta2)

    X1<-rep(X1,n2+1)
    X2<-rep(X2,each=n1+1)
    f1<- rep(f1,n2+1)
    f2<- rep(f2,each=n1+1)
    fall<- f1*f2

    pvals<-rep(NA,length(fall))
    for (i in 1:length(fall)){
        pvals[i]<-uncondExact2x2(X1[i],n1,X2[i],n2,...)$p.value
    }

    power<- sum( fall[ pvals<=alpha] )
    power
}