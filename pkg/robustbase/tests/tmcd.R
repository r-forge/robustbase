library(robustbase)

source(system.file("xtraR/test_MCD.R", package = "robustbase"))#-> doMCDdata
##          ../inst/xtraR/test_MCD.R
source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
## -> assertError(), relErr(), and:
showProc.time()

## -- now do it:
options(digits = 5)
set.seed(101) # <<-- sub-sampling algorithm now based on R's RNG and seed
doMCDdata()
doMCDdata(method="DetMCD"); warnings()
##                        vvvv no timing for 'R CMD Rdiff' outputs
doMCDdata(nrep = 12, time=FALSE)
doMCDdata(nrep = 12, time=FALSE, method="DetMCD"); warnings()
doMCDdata(nrep = 12, time=FALSE, method = "MASS")

###--- now the "close to singular" mahalanobis case:
set.seed(6)
(c3  <- covMcd(mort3))
(c3. <- covMcd(mort3, nsamp="deterministic"))
stopifnot(log(c3$crit) <= log(c3.$crit),
          print(log(c3.$crit / c3$crit)) <= 0.8)
## see 0.516 / 0.291 {with seed 7}
##
## rescale variables:
scaleV <- c(0.1, 0.1, 1, 1, .001, 0.1, 0.1, 100)
mm <- data.matrix(mort3) * rep(scaleV, each = nrow(mort3))
C3  <- covMcd(mm)
C3. <- covMcd(mm, nsamp="deterministic")
stopifnot(C3$mcd.wt == c3$mcd.wt)# here, not for all seeds!

## error ("computationally singular") with old (too high) default tolerance:
try( covMcd(mm, control= rrcov.control(tol = 1e-10)) )
try( covMcd(mm, control= rrcov.control(tol = 1e-10), nsamp="deterministic") )

showProc.time()

## "large" examples using different algo branches {seg.fault in version 0.4-4}:

n <- 600 ## - partitioning will be triggered
set.seed(1)
X <- matrix(round(100*rnorm(n * 3)), n, 3)
(cX  <- covMcd(X))
 cX. <- covMcd(X, nsamp="deterministic", scalefn = scaleTau2)
i <- names(cX); i <- i[!(i %in% c("call", "nsamp", "method", "raw.weights"))]
stopifnot(sum(cX.$raw.weights != cX$raw.weights) <= 2,
          all.equal(cX[i], cX.[i], tol= 1/9))

n <- 2000 ## - nesting will be triggered
set.seed(4)
X <- matrix(round(100*rnorm(n * 3)), n, 3)
set.seed(1)
summary(cX  <- covMcd(X)) # <- show newly activated  print.summary.mcd(.)
 cX. <- covMcd(X, nsamp="deterministic", scalefn = scaleTau2)
i2 <- i[i != "mcd.wt"]
stopifnot(print(sum(cX.$raw.weights != cX$raw.weights)) <= 3, # 2
          all.equal(cX[i2], cX.[i2], tol= 1/10))# 1/16

set.seed(1) ## testing of 'raw.only' :
cXo <- covMcd(X, raw.only=TRUE)
i <- paste0("raw.", c("cov", "center", "cnp2"))
stopifnot(cXo$raw.only, all.equal(cX[i], cXo[i], tol = 1e-15),
          c("best", "mah") %in% setdiff(names(cX), names(cXo)))
showProc.time()

## Now, some small sample cases:

## maximal values:
n. <- 10
p. <-  8
set.seed(44)
(X. <- cbind(1:n., round(10*rt(n.,3)), round(10*rt(n.,2)),
             matrix(round(10*rnorm(n. * (p.-3)), 1),  nrow = n., ncol = p.-3)))

## 2 x 1 ---> Error
r <- tryCatch(covMcd(X.[1:2, 2, drop=FALSE]), error=function(e)e)
stopifnot(inherits(r, "error"),
          grepl("too small sample size", r$message))

## 3 x 2 --- ditto
r <- tryCatch(covMcd(X.[1:3, 2:3]), error=function(e)e)
stopifnot(inherits(r, "error"),
          grepl("too small sample size", r$message))

## 5 x 3  [ n < 2 p  ! ]  --- also works for MASS
X <- X.[1:5, 1:3]
set.seed(101)
## the finite-sample correction is definitely doubtful:
summary(cc <- covMcd(X, use.correction = FALSE))
str(cc) ## best = 2 3 4 5
mcc <- MASS::cov.mcd(X)
stopifnot(cc$best == mcc$best,
          all.equal(cc$center, mcc$center, tolerance = 1e-10),
          all.equal(c(mcc$cov / cc$raw.cov), rep(0.673549282206, 3*3)))

## p = 4 -- 6 x 4 & 7 x 4  [ n < 2 p  ! ]
p <- 4
n <- 7
X <- X.[1:n, 1+(1:p)]
stopifnot(dim(X) == c(n,p))
(cc <- covMcd(X, use.correction = FALSE))
str(cc) ## best = 1 2 4 5 6 7
mcc <- MASS::cov.mcd(X)
stopifnot(cc$best == mcc$best,
          all.equal(cc$center, mcc$center, tolerance = 1e-10),
          all.equal(c(mcc$cov / cc$raw.cov), rep(0.7782486992881, p*p)))
n <- 6
X <- X[1:n,]
(cc <- covMcd(X, use.correction = FALSE))
mcc <- MASS::cov.mcd(X)
stopifnot(cc$best == mcc$best,
          all.equal(cc$center, mcc$center, tolerance = 1e-10),
          all.equal(c(mcc$cov / cc$raw.cov), rep(0.7528695976179, p*p)))

showProc.time()

## nsamp = "exact" -- here for p=7
coleman.x <- data.matrix(coleman[, 1:6])
showSys.time(CcX <- covMcd(coleman.x, nsamp= "exact"))
showSys.time(Ccd <- covMcd(coleman.x, nsamp= "deterministic"))
stopifnot(all.equal(CcX$best,
		    c(2, 5:9, 11,13, 14:16, 19:20), tolerance=0),
	  intersect(CcX$best, Ccd$best) == c(2,5,7,8,13,14,16,19,20),
          relErr(CcX$crit, Ccd$crit) < 0.35 # see ~ 0.34
)
summary(Ccd)

demo(determinMCD)
warnings()

