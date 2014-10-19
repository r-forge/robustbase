## VT::16.10.2014 - modified from tmcd.R - to test the
## deterministic option of covMcd.
##
library(robustbase)

source(system.file("xtraR/test_MCD.R", package = "robustbase"))#-> doMCDdata
##          ../inst/test_MCD.R

## -- now do it:
options(digits = 5)
set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed
doMCDdata(method="DETMCD")
##                vvvv no timing for 'R CMD Rdiff' outputs
doMCDdata(nrep = 12, time=FALSE, method="DETMCD")
doMCDdata(nrep = 12, time=FALSE, method = "MASS")

###--- now the "close to singular" mahalanobis case:
(c3 <- covMcd(mort3, nsamp="deterministic"))
## rescale variables:
scaleV <- c(0.1, 0.1, 1, 1, .001, 0.1, 0.1, 100)
mm <- data.matrix(mort3) * rep(scaleV, each = nrow(mort3))
C3 <- covMcd(mm, nsamp="deterministic")
stopifnot(C3$mcd.wt == c3$mcd.wt)
try(## error: with "old default tolerance:
  covMcd(mm, control= rrcov.control(tol = 1e-10), nsamp="deterministic")
)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''

## "large" examples using different algo branches {seg.fault in version 0.4-4}:
set.seed(1)

n <- 600 ## - partitioning will be triggered
X <- matrix(round(100*rnorm(n * 3)), n, 3)
cX <- covMcd(X, nsamp="deterministic")
cX
n <- 2000 ## - nesting will be triggered
X <- matrix(round(100*rnorm(n * 3)), n, 3)
cX <- covMcd(X, nsamp="deterministic")
cX

cat('Time elapsed: ', proc.time(),'\n')


## Now, some small sample cases:

## maximal values:
n. <- 10
p. <-  8
set.seed(44)
(X. <- cbind(1:n., round(10*rt(n.,3)), round(10*rt(n.,2)),
             matrix(round(10*rnorm(n. * (p.-3)), 1),  nrow = n., ncol = p.-3)))

## 2 x 1 ---> Error
r <- try(covMcd(X.[1:2, 2, drop=FALSE], nsamp="deterministic"), silent=TRUE)
stopifnot(inherits(r, "try-error"),
          grep("too small sample size", r) == 1)

## 3 x 2 --- ditto
r <- try(covMcd(X.[1:3, 2:3], nsamp="deterministic"), silent=TRUE)
stopifnot(inherits(r, "try-error"),
          grep("too small sample size", r) == 1)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
