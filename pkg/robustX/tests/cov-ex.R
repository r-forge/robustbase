library(robustX)
library(robustbase)
(newRB <- (packageVersion("robustbase") >= "0.99"))

sessionInfo()
packageDescription("robustX")
(ourBLAS <- grepl(print(normalizePath(R.home())),
                        normalizePath(extSoftVersion()[["BLAS"]]), fixed = TRUE))
## need extended precision (typically *includes* 64-bit):
doCheck <- (.Machine$sizeof.longdouble >= 16)
cat("doCheck (= have long double):", doCheck,"\n")

if(!dev.interactive(orNone=TRUE)) pdf("cov-ex.pdf")

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
    TRUE
}

UN <- function(L) lapply(L, unname)

chk.NN.new.old <- function(cNew, cNold, tol = 2e-15, tol.1 = 20*tol) {
    stopifnot(is.list(cNold$innc), length(n.i <- names(cNold$innc)) == 4)
    cat("classification accordance matrix:\n")
    print(table(new = cNew $classification,
                old = cNold$classification))
    report.stopifnot.all.eq(UN(cNew [1:4]),
                            UN(cNold[1:4]), tol=tol.1) &
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


try( # testing (tol=0 too small)
    chk.NN.new.old(cN, cN1, tol=0)
)
## This used to fail when we use R's instead of BLAS matrix products:
if(doCheck)
    chk.NN.new.old(cN, cN1, tol = 4e-15) # seen 1.1e-15 work


## for n = 500, you *do* see it
n <- 500
set.seed(12)
X <- rbwheel(n, 7, spherize=TRUE)

lattice::splom(X, cex=.1)
system.time(cNX1 <- covNN.1(X))# 0.82  0.273
system.time(cNX  <- covNNC (X))# 0.66  0.163
system.time(cM   <- covMcd (X))# 0.151 0.097  <- !
# NB: *slower* times above, when using R's instead of BLAS matrix prod

try( # --> show relative difference(s):
    chk.NN.new.old(cNX, cNX1, tol=0)
)
if(doCheck && ourBLAS) # did fail with ATLAS in R-devel 2023-1-1
    chk.NN.new.old(cNX, cNX1)

stopifnot(exprs = {
    all.equal(1900.4208,   kappa(cM $cov))# 1990.8.. then  1900.421
    all.equal(4.485807117, kappa(cNX$cov))
    all.equal(1.047781251, kappa(cov(X)))
})

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

(C.mcd <- covMcd(X1)$cov)
Cm <- as.matrix(if(newRB) 4.8848 else 7.790004)
all.equal(Cm, C.mcd, tol=0) # 6.633e-6
stopifnot(all.equal(Cm, C.mcd, tol = 2e-5))


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


## Less trivial data  --- also used in ../man/BACON.Rd
data(starsCYG, package = "robustbase")

op <- options(warn = 2)# no warnings allowed
str(B.ST <- with(starsCYG, BACON(x = log.Te, y = log.light)))
(Bgood <- which(B.ST$subset))
.Platform$r_arch
## 32-bit <-> 64-bit different results {tested on Linux & Windows Server}
is32 <- .Machine$sizeof.pointer == 4 ## <- should work for Linux/MacOS/Windows
isMac <- Sys.info()[["sysname"]] == "Darwin"
isSun <- Sys.info()[["sysname"]] == "SunOS"
isWin <- .Platform$OS.type == "windows"
(Platf_arch <- paste0(.Platform$r_arch, # maybe distinguish more, e.g. "no-ldouble", "BLAS" version ,  ??
                      if(isWin) "_win"))
knownPl <- Platf_arch %in% c("x64", "i386_win", "i386") # update!
stopifnot(exprs = {
    switch(
        Platf_arch,
        "i386_win" =, # Platform: i386-w64-mingw32/i386 (32-bit); Windows Server x64 (build 14393)
        "x64"  = identical(Bgood, c(25L, 27:29, 33L,        38L,      43L, 45L)),
        "i386" = identical(Bgood, c(25L, 27:29, 33L, 35:36, 38L, 42L, 43L, 45L)),
                                        # older version of i386-windows
        ## other platforms:
        {
            message("Platform architecture (see above) not yet tested for BACON result")
            TRUE
        })
})


plot(starsCYG)
lmST <- lm(log.light ~ log.Te, data = starsCYG)
abline(lmST, col = "gray") # least squares line
## 'subset': A good set of of points (to determine regression):
colB <- adjustcolor(2, 1/2)
points(log.light ~ log.Te, data = starsCYG, subset = B.ST$subset,
       pch = 19, cex = 1.5, col = colB)
## A BACON-derived line:
lmB <- lm(log.light ~ log.Te, data = starsCYG, subset = B.ST$subset)
abline(lmB, col = colB, lwd = 2)
cfT <- switch(Platf_arch # use lm(..., subset = sfsmisc::inverseWhich(Bgood))
            , "x64"  = c(-10.130429, 3.4450151)
            , "i386" = c(-9.5759001, 3.3109007)
             ## otherwise:
            , NULL)
cf <- unname(coef(lmB))
dput(signif(cf, 8))
if(!is.null(cfT)) withAutoprint({
    all.equal(cf, cfT, tol=0)# 64b: 6.4e-10
    stopifnot(all.equal(cf, cfT))
})

options(op)
