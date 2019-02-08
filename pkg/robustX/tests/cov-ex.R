library(robustX)
library(robustbase)

covNN.1 <- robustX:::covNNC1  ## the original definition (2003)

data(iris)
system.time(cN1 <- covNN.1(iris[-5]))
system.time(cN  <- covNNC (iris[-5]))# faster indeed

## report.and.stop.if.not.all.equal
report.stopifnot.all.eq <- function(a,b, tol, ...) {
    call <- sys.call()
    ae <- all.equal(a,b, tol=tol, ...)
    call[[1]] <- quote(all.equal)
    if(!isTRUE(ae))
	stop(sprintf("Not %s:\n%s\n\n", deparse(call),
		     paste(ae, collapse="\n")),
	     call.=FALSE)
    ## else
    invisible(TRUE)
}

UN <- function(L) lapply(L, unname)

chk.NN.new.old <- function(cNew, cNold, tol = 2e-15, tol.1 = 20*tol) {
    stopifnot(is.list(cNold$innc), length(n.i <- names(cNold$innc)) == 4)
    cat("classification accordance matrix:\n")
    print(table(new = cNew $classification,
                old = cNold$classification))
    report.stopifnot.all.eq(UN(cNew [1:4]),
                            UN(cNold[1:4]), tol=tol.1)
    report.stopifnot.all.eq(cNew $innc[n.i],
                            cNold$innc[n.i], tol=tol)
}

summ.NN <- function(cNN, digits = 3) {
    cbind(class = cNN$classification,
          pprob = round(cNN$postprob, digits),
          incc.p= round(cNN$innc$postprob, digits))
}

s1 <- summ.NN(cN1)
ss <- summ.NN(cN)
if(isTRUE(all.equal(ss, s1))) ss else cbind(ss, s1)


try( # testing
    chk.NN.new.old(cN, cN1, tol=0)
)
## need extended precision (typically *includes* 64-bit):
doCheck <- (.Machine$sizeof.longdouble >= 16)
cat("doCheck (= have long double):", doCheck,"\n")

## This fails (interestingly) when we use R's instead of BLAS matrix products:
'MM:  no it now works ! ??'
if(doCheck) try( chk.NN.new.old(cN, cN1) )



## for n = 500, you *do* see it
n <- 500
set.seed(12)
X <- rbwheel(n, 7, spherize=TRUE)

lattice::splom(X, cex=.1)
system.time(cNX1 <- covNN.1(X))# 0.82  0.273
system.time(cNX  <- covNNC (X))# 0.66  0.163
system.time(cM   <- covMcd (X))# 0.151 0.097  <- !
# NB: *slower* times above, when using R's instead of BLAS matrix prod

try( # testing
    chk.NN.new.old(cNX, cNX1, tol=0)
)
if(doCheck)
    chk.NN.new.old(cNX, cNX1)

kappa(cM $cov)# 1990.8.. then  1900.421
kappa(cNX$cov)#    4.4858
kappa(cov(X)) #    1.0478

## ---- d = 1 :
X1 <- cbind(c(1:6, 1000))

var(X1)
##          [,1]
## [1,] 141861.8
## if 1000 was not an outlier:
var(1:7) ## 4.666667

covNNC(X1)$cov ## -- really not at all robust:
##          [,1]
## [1,] 121595.8

C.mcd <- covMcd(X1)$cov
##          [,1]
## [1,] 7.790004
all.equal(C.mcd, as.matrix(7.79), tol=0)
stopifnot(all.equal(C.mcd, as.matrix(7.79), tol = 1e-6))


MASS::cov.rob(X1)$cov
##      [,1]
## [1,]  3.5
(C.B <- BACON(X1)$cov)
##      [,1]
## [1,]  3.5
all.equal(C.B, as.matrix(3.5), tol=0)
stopifnot(all.equal(C.B, as.matrix(3.5)))

if(FALSE) ## FIXME (in robustbase!): should work for  p=1
    covOGK(X1)$cov
