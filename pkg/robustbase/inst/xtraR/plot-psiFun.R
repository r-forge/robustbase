## Functions to plot and check psi-functions
## used in ../../tests/lmrob-psifns.R,
##	   ../../tests/psi-rho-etc.R
##     and ../doc/psi_functions.Rnw  vignette

psiF <- robustbase:::lmrob.psifun # deriv = -1 (rho), 0, 1
chiF <- robustbase:::lmrob.chifun # rho(.) normalized to max|.| = 1;  deriv
wgtF <- robustbase:::lmrob.wgtfun

## Original Author of functions: Martin Maechler, Date: 13 Aug 2010, 10:17
plot.psiFun <- function(x, m.psi, psi, par, shortMain = FALSE,
                        col = c("black", "red3", "blue3", "dark green"),
                        leg.loc = "right", ...) {
    fExprs <- if (ncol(m.psi) > 4) {
        quote(list(rho(x), psi(x), {psi*minute}(x), w(x) == psi(x)/x,
               {w*minute}(x)))
    } else quote(list(rho(x), psi(x), {psi*minute}(x), w(x) == psi(x)/x))
    tit <- if(shortMain)
	substitute(rho(x) ~ "etc, with" ~ psi*"-type" == PSI(PPP),
		   list(PSI = psi, PPP = paste(formatC(par), collapse=",")))
    else
	substitute(FFF ~~ ~~ " with "~~ psi*"-type" == PSI(PPP),
		   list(FFF = fExprs, PSI = psi,
			PPP = paste(formatC(par), collapse=",")))
    matplot(x, m.psi, col=col, lty=1, type="l", main = tit,
            ylab = quote(f(x)), xlab = quote(x), ...)
    abline(h=0,v=0, lty=3, col="gray30")
    fE <- fExprs; fE[[1]] <- as.name("expression")
    legend(leg.loc, inset=.02, eval(fE), col=col, lty=1, bty="n")
    invisible(cbind(x=x, m.psi))
}

p.psiFun <- function(x, psi, par, ...)
{
    m.psi <- cbind(rho    = psiF(x, par, psi,deriv=-1),
                   psi    = psiF(x, par, psi,deriv= 0),
                   "psi'" = psiF(x, par, psi,deriv= 1),
                   wgt    = wgtF(x, par, psi))
    plot.psiFun(x, m.psi, psi, par, ...)
}
p.psiFun2 <- function(x, psi, par, ...)
    p.psiFun(x, psi, par, shortMain = TRUE,
             leg.loc = "bottomright", ylim = c(-3, 5))
## the same for objects of psi_func-class
p.psiFun3 <- function(x, object, ...)
{
    ## Author: Martin Maechler, Date: 13 Aug 2010, 10:17
    m.psi <- cbind(rho    = object@rho(x),
                   psi    = object@psi(x),
                   "psi'" = object@Dpsi(x),
                   wgt    = object@wgt(x),
                   "wgt'" = object@Dwgt(x))
    plot.psiFun(x, m.psi, object@name, unlist(formals(object@rho)[-1]),
                col = c("black", "red3", "blue3", "dark green", "light green"),
                ...)
}


mids <- function(x) (x[-1]+x[-length(x)])/2
chkPsiDeriv <- function(m.psi, tol = 1e-4) {
    ## m.psi: matrix as from p.psiFun()
    stopifnot(length(tol) > 0, tol >= 0,
              is.numeric(psi <- m.psi[,"psi"]),
              is.numeric(dx  <- diff(x <- m.psi[,"x"])))
    if(length(tol) < 2) tol[2] <- 10*tol[1]
    xn0 <- abs(x) > 1e-5
    c(all.equal(mids(psi), diff(m.psi[,"rho"])/dx, tol=tol[1]), # rho'  == psi
      all.equal(mids(m.psi[,"psi'"]), diff(psi)/dx, tol=tol[2]),# psi'  == psip
      all.equal(m.psi[xn0,"wgt"], (psi/x)[xn0], tol= tol[1]/10))# psi/x == wgt
}
