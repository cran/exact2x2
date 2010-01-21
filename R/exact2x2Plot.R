exact2x2Plot<-function(x,
    OR=NULL,
    ndiv=1000,
    method="minlike",
    orRange=c(.01,10),
    dolog=TRUE,
    dolines=FALSE,
    dopoints=FALSE,
    newplot=TRUE,...){

    #m<-sum(x[1,])
    #n<-sum(x[2,])
    #k<-sum(x[,1])
    #x<-x[1,1]

    if (dolog & is.null(OR)){
        minor<-log10(orRange[1])
        maxor<-log10(orRange[2])
        OR<-10^(minor + (0:ndiv)*(maxor-minor)/ndiv )
    } else if (is.null(OR)){
        minor<-(orRange[1])
        maxor<-(orRange[2])
        OR<-(minor + (0:ndiv)*(maxor-minor)/ndiv )
    }
    p<-exact2x2Pvals(x,OR,method=method)$pvals

    if (newplot){
        if (dolog){ LOG<-"x"
        } else LOG<-""
        plot(OR,p,log=LOG,xlab="Null Hypothesis Odds Ratio",ylab="two-sided p-value",...)
    } 
    if (dopoints){
        points(OR,p,...)
    } else if (dolines){
        lines(OR,p,...)
    }
}
