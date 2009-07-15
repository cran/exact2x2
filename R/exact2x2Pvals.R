`exact2x2Pvals` <-
function(x,or,relErr=1+10^(-7),method="minlike"){

    m<-sum(x[1,])
    n<-sum(x[2,])
    k<-sum(x[,1])
    x<-x[1,1]



    if (method=="fisher"){
        method<-"minlike"
        warning("method='fisher' changed to 'minlike' ")
    }
    nor<-length(or)
    #fisher<-blaker<-central<-rep(NA,nor)
    pvals<-rep(NA,nor)
    support<- max(0,k-n):min(m,k)
    X<- support==x
    ns<-length(support)
    logf0 <- dhyper(support, m, n, k, log = TRUE)
    ## use dnhyper and relErr as in fisher.test
    dnhyper <- function(OR) {
        d <- logf0 + log(OR) * support
        d <- exp(d - max(d))
        d/sum(d)
    }
    if (method=="blaker"){
        for (i in 1:nor){
            f<-dnhyper(or[i])
            d<-f[X]
            #fisher[i]<-sum(f[f<=d*relErr])
            F<-cumsum(f)
            Fbar<-cumsum(f[ns:1])[ns:1]
            lower<-F[X]
            upper<-Fbar[X]
            if (lower<=upper){
                if (any(Fbar<=lower*relErr)){
                    pvals[i]<-lower + max(Fbar[Fbar<=lower*relErr])
                } else pvals[i]<-lower
            } else if (upper<lower){
                if (any(F<=upper*relErr)){
                    pvals[i]<-upper + max(F[F<=upper*relErr])
                } else pvals[i]<-upper
            } 
        }
    } else if (method=="minlike"){
        for (i in 1:nor){
            f<-dnhyper(or[i])
            d<-f[X]
            pvals[i]<-sum(f[f<=d*relErr])
        }   
    } 
    if (method=="all" | method=="central"){
        Xlo<-support<=x
        Xhi<-support>=x
        for (i in 1:nor){
            f<-dnhyper(or[i])
            pvals[i]<-2*min(sum(f[Xlo]),sum(f[Xhi]))
        }
        pvals<-pmin(1,pvals)
    }
    out<-list(or=or,pvals=pvals)
    out
}

