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
    ret <- lmrob.split(mf)
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
m_s_subsample <- function(x1, x2, y, control, orthogonalize=TRUE) {
    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)
    y <- y
    storage.mode(x1) <- "double"
    storage.mode(x2) <- "double"
    storage.mode(y) <- "double"
    
    z <- .C(robustbase:::R_lmrob_M_S,
            X1=x1,
            X2=x2,
            y=y,
            n=length(y),
            p1=ncol(x1),
            p2=ncol(x2),
            nResample=as.integer(control$nResample),
            scale=double(1),
            b1=double(ncol(x1)),
            b2=double(ncol(x2)),
            tuning_chi=as.double(control$tuning.chi),
            ipsi=as.integer(lmrob.psi2ipsi(control$psi)),
            bb=as.double(control$bb),
            K_m_s=as.integer(control$k.m_s),
            max_k=as.integer(control$max.k),
            rel_tol=as.double(control$rel.tol),
            converged=logical(1),
            trace_lev=as.integer(control$trace.lev),
            do_descent=FALSE,
            orthogonalize=as.logical(orthogonalize))
    z[c("b1", "b2", "scale")]
}

control <- lmrob.control()
obj <- lm(Y ~ Region + X1 + X2 + X3, education)
splt <- lmrob.split(obj$model)
y <- education$Y

## test orthogonalizing
x1 <- splt$x1
x2 <- splt$x2
tmp <- lmrob.lar(x1, y, control$rel.tol)
y.tilde <- tmp$resid
t1 <- tmp$coef
x2.tilde <- x2
T2 <- matrix(0, nrow=ncol(x1), ncol=ncol(x2))
for (i in 1:ncol(x2)) {
    tmp <- lmrob.lar(x1, x2[,i], control$rel.tol)
    x2.tilde[,i] <- tmp$resid
    T2[,i] <- tmp$coef
}
set.seed(10)
res1 <- m_s_subsample(x1, x2.tilde, y.tilde, control, FALSE)
res1 <- within(res1, b1 <- drop(t1 + b1 - T2 %*% b2))
set.seed(10)
res2 <- m_s_subsample(x1, x2, y, control, TRUE)
stopifnot(all.equal(res1, res2))

res <- list()
set.seed(0)
for (i in 1:100) {
    tmp <- m_s_subsample(x1, x2.tilde, y.tilde, control, FALSE)
    res[[i]] <- unlist(within(tmp, b1 <- drop(t1 + b1 - T2 %*% b2)))
}
res <- do.call(rbind, res)
## show a summary of the results
summary(res)
## compare with fast S solution
obj <- lmrob(Y ~ Region + X1 + X2 + X3, education, init="S")
coef(obj)
obj$scale
