fisherCalcMidp<-function(x,or,alternative,conf.int,conf.level){
    # much of the code is copied from fisher.test
    # we just need to make some modifications for the midp
    # adjustment
    m <- sum(x[, 1L])
    n <- sum(x[, 2L])
    k <- sum(x[1L, ])
    x <- x[1L, 1L]
    lo <- max(0L, k - n)
    hi <- min(k, m)
    support <- lo:hi
    logdc <- dhyper(support, m, n, k, log = TRUE)
    dnhyper <- function(ncp) {
        d <- logdc + log(ncp) * support
        d <- exp(d - max(d))
        d/sum(d)
    }
    # modify pnhyper from fisher.test to pnhyperMidp
    # usually just subtract off 0.5*dnhyper(or)[support==x] 
    # 
    pnhyperMidp <- function(q, ncp = 1, upper.tail = FALSE) {
        if (ncp == 1) {
            return(if (upper.tail) phyper(x - 1, m, n, k, 
              lower.tail = FALSE) - 0.5*dhyper(x,m,n,k) else phyper(x, 
              m, n, k)- 0.5*dhyper(x,m,n,k))
        }
        if (ncp == 0) {
            return(as.numeric(if (upper.tail) q <= lo else q >= 
              lo) - 0.5*as.numeric(q==lo))
        }
        if (ncp == Inf) {
            return(as.numeric(if (upper.tail) q <= hi else q >= 
              hi) - 0.5*as.numeric(q==hi))
        }
        sum(dnhyper(ncp)[if (upper.tail) support >= q else support <= 
                q]) - 0.5*dnhyper(ncp)[support==q]
    }

    pless<-  pnhyperMidp(x, or)
    pgreater<- pnhyperMidp(x, or, upper.tail = TRUE)
    # only do mid-p for central two-sided p-values
    PVAL <- switch(alternative, 
        less =pless, 
        greater = pgreater, 
        two.sided = min(1,2*min(pless,pgreater)))
    if (conf.int) {
        ncp.U <- function(x, alpha) {
            if (x == hi) return(Inf)
            if (x==lo & alpha>=0.5) return(0)
            p <- pnhyperMidp(x, 1)
            if (p < alpha) 
              uniroot(function(t) pnhyperMidp(x, t) - alpha, 
                c(0, 1))$root
            else if (p > alpha) 
              1/uniroot(function(t) pnhyperMidp(x, 1/t) - alpha, 
                c(.Machine$double.eps, 1))$root
            else 1
        }
        ncp.L <- function(x, alpha) {
           if (x == lo) return(0)
           if (x==hi & alpha>=0.5) return(Inf)
           p <- pnhyperMidp(x, 1, upper.tail = TRUE)
            if (p > alpha) 
              uniroot(function(t) pnhyperMidp(x, t, upper.tail = TRUE) - 
                alpha, c(0, 1))$root
            else if (p < alpha) 
              1/uniroot(function(t) pnhyperMidp(x, 1/t, upper.tail = TRUE) - 
                alpha, c(.Machine$double.eps, 1))$root
            else 1
        }
        CINT <- switch(alternative, 
            less = c(0, ncp.U(x, 1 - conf.level)), 
            greater = c(ncp.L(x, 1 - conf.level), Inf), 
            two.sided = {
                alpha <- (1 - conf.level)/2
                c(ncp.L(x, alpha), ncp.U(x, alpha))
        })
        attr(CINT, "conf.level") <- conf.level
    }
    out<-  list(conf.int = if (conf.int) CINT, p.value=PVAL)
    out
}


#check
#x<-matrix(c(4,0,12,5),2,2)
#x<-matrix(c(0,4,12,5),2,2)
#x<-matrix(c(2,2,12,5),2,2)
#r<-fisherCalcMidp(x,1,"less",TRUE,.95)
#r2<-fisherCalcMidp(x,1,"less",TRUE,1-r$p.value)
#r2

