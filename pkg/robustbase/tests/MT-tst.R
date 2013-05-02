require("robustbase")

##---> ./poisson-ex.R
##     ~~~~~~~~~~~~~~  for more glmrobMT() tests

source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
## -> assertError(), showSys.time(), ...
## Newer version of the above test-tools-1.R contain this:
assert.EQ <- function(target, current, tol = if(show) 0 else 1e-15,
                      show = FALSE, ...) {
    ## Purpose: check equality *and* show non-equality
    ## ----------------------------------------------------------------------
    ## show: if TRUE, return (and hence typically print) all.equal(...)
    if(show) all.equal(target, current, tol = tol)
    else if(!isTRUE(r <- all.equal(target, current, tol = tol)))
	stop("all.equal() |-> ", paste(r, collapse=sprintf("%-19s","\n")))
}


(doExtras <- robustbase:::doExtras())

## Explore the espRho() function: ---------------------------------------------
pdf("MT-E_rho.pdf")
E.rho <- robustbase:::espRho
lambdas <- ((1:10)/2)^2
cws <- c(1, 1.5, 1.75, 2, 2.25, 3)
(gr <- expand.grid(lam = lambdas, cw = cws))

Egr <- apply(gr, 1, function(r) {
    lam <- r[["lam"]]; cw <- r[["cw"]]; sL <- sqrt(lam)
    xx <- seq(lam - 2*sL, lam + 2*sL, length=17)
    vapply(xx, function(X) E.rho(lam, xx=X, cw=cw), NA_real_)
})
str(Egr)# 17 x 60
mLeg <- function(pos, type="o")
    legend(pos, legend=paste("lambda = ", format(lambdas, digits=2)),
           lty=1:5, col=1:6, pch= c(1:9, 0, letters, LETTERS), bty="n")
matplot(Egr[, gr[,"cw"] == 1.0 ], type="o", main= "c_w = 1.0" ); mLeg("bottomright")
matplot(Egr[, gr[,"cw"] == 1.5 ], type="o", main= "c_w = 1.5" ); mLeg("bottomright")
matplot(Egr[, gr[,"cw"] == 1.75], type="o", main= "c_w = 1.75"); mLeg("bottomright")
matplot(Egr[, gr[,"cw"] == 2.0 ], type="o", main= "c_w = 2.0" ); mLeg("bottomright")
matplot(Egr[, gr[,"cw"] == 2.25], type="o", main= "c_w = 2.25"); mLeg("bottomright")
matplot(Egr[, gr[,"cw"] == 3.0 ], type="o", main= "c_w = 3.0" ); mLeg("bottomright")

dev.off()


## Explore the m() function: ---------------------------------------------
pdf("MT-m_rho.pdf")

mkM <- robustbase:::mk.m_rho
m21 <- mkM(2.1, recompute=TRUE)# the default 'cw = 2.1'
m16 <- mkM(1.6, recompute=TRUE)
p.m2 <- function(mrho, from = 0, to, col=2, addKnots=TRUE, pchK=4, cexK=1.5, ...) {
    stopifnot(is.function(mrho))
    curve(mrho, from, to, col=col, ...)
    curve(sqrt(x), add=TRUE, col=adjustcolor("gray",.5), lwd=2)
    if(addKnots) { e <- environment(mrho); points(e$x0, e$y0, pch=pchK, cex=cexK) }
}
p.m.diff <- function(mrho, from = 0, to, col=2, addKnots=TRUE, pchK=4, cexK=1.5, ...) {
    stopifnot(is.function(mrho))
    curve(mrho(x) - sqrt(x), from=from, to=to, n=512, col=col, ...)
    abline(h=0,lty=3)
    if(addKnots) {
        e <- environment(mrho); x <- e$x0
        if(is.numeric(x))
            points(x, e$y0 - sqrt(x), pch=pchK, cex=cexK)
        else warning("'addKnots' not available: No knots in function's environment")
    }
}

p.m2(m21, to=10)
p.m2(m16, to=10)
p.m2(m21, to=50)
p.m2(m21, to=120, cexK=.8)
p.m.diff(m21, to=120, cex=.5)# pchK="."
p.m.diff(m16, to=120, cex=.5)# pchK="."

mm21 <- function(.) robustbase:::mm(., m21)
environment(mm21) <- environment(m21)# <- for p.m()
p.m2(mm21, to=120, cexK=.8)
p.m.diff(mm21, to=120, cexK=.8)#-- discontinuity at 100 !!
## TODO: ways to improve!

## TODO first: look at more cases (cw)

dev.off()
##-------------------------------------------------------- end m(.) -------------


## The simple intercept example from  ./glmrob-1.R
set.seed(113)
y <- rpois(17, lambda = 4)
y[1:2] <- 99:100 # outliers
y.1 <- y
## To call Victor's version of glmrobMT()  for this case, we need  "hoop jumps"
x.1 <- matrix(0, nrow=length(y.1), ncol=0)

options("robustbase:m_rho_recompute" = TRUE)#-> recompute in any case:
showSys.time( r <- glmrobMT(x.1, y.1, nsubm=100) )# some output
str(r)

## was   c(ini = 1.30833281965018, est = 1.29369680430613)
## then  c(ini = 1.30833281965018, est = 1.29369680422799)
##       c(ini = 1.30833281965018, est = 1.29369680430627)
r.64b <- c(ini = 1.30833281965018, est = 1.29369680452016)
stopifnot(r$converged)
assert.EQ(r$initial, r.64b[["ini"]], tol = 1e-13)# rel.diff: 3.394.e-16
assert.EQ(r$final,   r.64b[["est"]], tol = 1e-09)# as long we use different optim())


## now, as the algorithm has a random start:
set.seed(7)
nSim <- if(doExtras) 20 else 2
showSys.time(LL <- replicate(nSim, glmrobMT(x.1, y.1, trace.lev=0),
                             simplify=FALSE))
ini <- sapply(LL, `[[`, "initial")
est <- sapply(LL, `[[`, "final")
## surprise:  all the 20 initial estimators are identical:
stopifnot(diff(range(ini)) == 0,
          diff(range(est)) == 0)
## probably too accurate ... but ok, for now
assert.EQ(est[1], r.64b[["est"]], tol = 1e-13)
assert.EQ(ini[1], r.64b[["ini"]], tol = 1e-13)



cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
## "Platform" info
SysI <- Sys.info()[c(1:2,4:5)]
if(require("sfsmisc")) c(SysI, MIPS=Sys.MIPS(), Sys.sizes()) else SysI
