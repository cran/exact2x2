`exact2x2` <-
function(x, y = NULL, or = 1, alternative = "two.sided", 
    tsmethod="minlike", 
    conf.int = TRUE, conf.level = 0.95, 
    tol=0.00001,conditional=TRUE) 
{
    if (!conditional) stop("unconditional exact tests not supported")
    if (tol<.Machine$double.eps^.5) warning("tol set very small, make sure pnhyper in exact2x2CI is this accurate")
    # copied setup code from fisher.test, july 1,2009
    DNAME <- deparse(substitute(x))
    METHOD <- "Fisher's Exact Test for Count Data"
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (is.matrix(x)) {
        if (any(dim(x) < 2)) 
            stop("'x' must have at least 2 rows and columns")
        if (!is.numeric(x) || any(x < 0) || any(is.na(x))) 
            stop("all entries of 'x' must be nonnegative and finite")
        if (!is.integer(x)) {
            xo <- x
            x <- round(x)
            if (any(x > .Machine$integer.max)) 
                stop("'x' has entries too large to be integer")
            if (!identical(TRUE, (ax <- all.equal(xo, x)))) 
                warning("'x' has been rounded to integer: ", 
                  ax)
            storage.mode(x) <- "integer"
        }
    }
    else {
        if (is.null(y)) 
            stop("if 'x' is not a matrix, 'y' must be given")
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        DNAME <- paste(DNAME, "and", deparse(substitute(y)))
        OK <- complete.cases(x, y)
        x <- factor(x[OK])
        y <- factor(y[OK])
        if ((nlevels(x) < 2) || (nlevels(y) < 2)) 
            stop("'x' and 'y' must have at least 2 levels")
        x <- table(x, y)
    }
    alternative <- char.expand(alternative, c("two.sided", 
        "less", "greater"))
    if (length(alternative) > 1 || is.na(alternative)) 
        stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
    # end of setup code from fisher.test
    nr <- nrow(x)
    nc <- ncol(x)
    if (nr!=2 || nc!=2) stop("table must be 2 by 2")
    if (alternative=="less" | alternative=="greater"){
        OUT<-fisher.test(x,or=or,alternative=alternative,conf.int=conf.int,conf.level=conf.level)
        OUT$method<-"One-sided Fisher's Exact Test"
    } else {
        xmat<-x
        m <- sum(x[, 1])
        n <- sum(x[, 2])
        k <- sum(x[1, ])
        x <- x[1, 1]
        tsmethod<-char.expand(tsmethod,c("minlike","blaker","central"))
        if (tsmethod!="blaker" & tsmethod!="minlike" & tsmethod!="central") stop("tsmethod must be 'blaker', 'central' or 'minlike'.")
        if (tsmethod=="minlike"){
            OUT<-fisher.test(xmat,or=or,alternative="two.sided",conf.int=FALSE) 
            if (conf.int){
                OUT$conf.int<-exact2x2CI(xmat,method="minlike",conf.level=conf.level,
                    tol=tol)
            } else {
                OUT$conf.int<-NULL
            }
            OUT$method<-"Two-sided Fisher's Exact Test (usual method using minimum likelihood)"
        } else if (tsmethod=="blaker"){
            OUT<-fisher.test(xmat,or=or,alternative="two.sided",conf.int=FALSE) 
            if (conf.int){
                OUT$conf.int<-exact2x2CI(xmat,method="blaker",conf.level=conf.level,
                    tol=tol)
            } else {
                OUT$conf.int<-NULL
            }
            OUT$p.value<-exact2x2Pvals(xmat,or,method="blaker")$pvals
            OUT$method<-"Blaker's Exact Test"
        } else if (tsmethod=="central"){
            OUT<-fisher.test(xmat,or=or,alternative=alternative,conf.int=conf.int,
                conf.level=conf.level)
            pless<-fisher.test(xmat,or=or,alternative="less",conf.int=FALSE)$p.value
            pgreater<-fisher.test(xmat,or=or,alternative="greater",conf.int=FALSE)$p.value
            OUT$p.value<-min(1,2*min(pless,pgreater))
            OUT$method<-"Central Fisher's Exact Test"
        } 
    }
    OUT$data.name<-DNAME
    OUT
}