## Test implementation of M-S estimator
require(robustbase)
lmrob.conv.cc <- robustbase:::lmrob.conv.cc
lmrob.psi2ipsi <- robustbase:::lmrob.psi2ipsi

## dataset with factors and continuous variables:
data(education)
education <- within(education, Region <- factor(Region))
## for testing purposes:
education2 <- within(education, Group <- factor(rep(1:3, length.out=length(Region))))
                     
## Test lmrob.split
testFun <- function(formula, x1.idx) {
    obj <- lm(formula, education2)
    mf <- obj$model
    mt <- attr(mf, "terms")
    x <- model.matrix(mt, mf) ## , contrasts)
    y <- model.response(mf, "numeric")
    ret <- lmrob.split(x, y, lmrob.control(), mf)
    if (missing(x1.idx)) {
        print(ret$x1.idx)
        return(which(unname(ret$x1.idx)))
    }
    stopifnot(all.equal(x1.idx, which(unname(ret$x1.idx))))
}
testFun(Y ~ 1, integer(0))
testFun(Y ~ X1*X2*X3, integer(0))
testFun(Y ~ Region + X1 + X2 + X3, 1:4)
testFun(Y ~ 0 + Region + X1 + X2 + X3, 1:4)
testFun(Y ~ Region*X1 + X2 + X3, c(1:5, 8:10))
testFun(Y ~ Region*X1 + X2 + X3 + Region*Group, c(1:5, 8:18))
testFun(Y ~ Region*X1 + X2 + X3 + Region*Group*X2, c(1:6, 8:29))
testFun(Y ~ Region*X1 + X2 + Region*Group*X2, 1:28)
testFun(Y ~ Region*X1 + X2 + Region:Group:X2, 1:21)
testFun(Y ~ Region*X1 + X2*X3 + Region:Group:X2, c(1:6, 8:10, 12:23))
testFun(Y ~ (X1+X2+X3+Region)^2, c(1:7,10:12,14:19))
testFun(Y ~ (X1+X2+X3+Region)^3, c(1:19, 21:29))
testFun(Y ~ (X1+X2+X3+Region)^4, 1:32)

## Test subsampling algorithm
m_s_subsample <- function(splt, control) {
    x1 <- as.matrix(splt$x1)
    x2.tilde <- as.matrix(splt$x2.tilde)
    y.tilde <- splt$y.tilde
    storage.mode(x1) <- "double"
    storage.mode(x2.tilde) <- "double"
    storage.mode(y.tilde) <- "double"
    
    z <- .C(robustbase:::R_lmrob_M_S,
            X1=x1,
            X2=x2.tilde,
            y=y.tilde,
            n=length(y.tilde),
            p1=ncol(x1),
            p2=ncol(x2.tilde),
            nResample=as.integer(control$nResample),
            scale=double(1),
            b1=double(ncol(x1)),
            b2=double(ncol(x2.tilde)),
            tuning_chi=as.double(control$tuning.chi),
            ipsi=as.integer(lmrob.psi2ipsi(control$psi)),
            bb=as.double(control$bb),
            K_m_s=as.integer(control$k.m_s),
            max_k=as.integer(control$max.k),
            rel_tol=as.double(control$rel.tol),
            converged=logical(1),
            trace_lev=as.integer(control$trace.lev),
            do.descent=FALSE)

    b1 <- drop(splt$t1 + z$b1 - splt$T2 %*% z$b2)
    list(b1=b1, b2=z$b2, scale=z$scale)
}

control <- lmrob.control()
obj <- lm(Y ~ Region + X1 + X2 + X3, education)
splt <- lmrob.split(model.matrix(obj), education$Y, control, obj$model)
res <- list()
set.seed(0)
for (i in 1:100)
    res[[i]] <- unlist(m_s_subsample(splt, control))
res <- do.call(rbind, res)
## show a summary of the results
summary(res)
## compare with fast S solution
obj <- lmrob(Y ~ Region + X1 + X2 + X3, education, init="S")
coef(obj)
obj$scale
