## Functions to plot and check psi-functions
## used in ../../tests/lmrob-psifns.R,
##	   ../../tests/psi-rho-etc.R
##     and ../doc/psi_functions.Rnw  vignette


## Original Author of functions: Martin Maechler, Date: 13 Aug 2010, 10:17

p.psiFun <- function(x, psi, par, main=FALSE, ...)
{
    m.psi <- cbind(rho  = .M.psi(x, par, psi,deriv=-1),
                   psi  = .M.psi(x, par, psi,deriv= 0),
                   Dpsi = .M.psi(x, par, psi,deriv= 1),
                   wgt  = .M.wgt(x, par, psi))
    robustbase:::matplotPsi(x, m.psi, psi=psi, par=par, main=main, ...)
}
p.psiFun2 <- function(x, psi, par, main="short", ...)
    p.psiFun(x, psi, par, main=main, leg.loc= "bottomright", ylim = c(-2.2, 6))
## for psi_func class objects: simply use plot() method.

mids <- function(x) (x[-1]+x[-length(x)])/2

##' @title Check consistency of psi/chi/wgt/.. functions
##' @param m.psi matrix as from p.psiFun()
##' @param tol
##' @return concatenation of \code{\link{all.equal}} results
##' @author Martin Maechler
chkPsiDeriv <- function(m.psi, tol = 1e-4) {
    stopifnot(length(tol) > 0, tol >= 0,
              is.numeric(psi <- m.psi[,"psi"]),
              is.numeric(dx  <- diff(x <- m.psi[,"x"])))
    if(length(tol) < 2) tol[2] <- 10*tol[1]
    xn0 <- abs(x) > 1e-5
    c(all.equal(mids(psi), diff(m.psi[,"rho"])/dx, tol=tol[1]), # rho'  == psi
      all.equal(mids(m.psi[,"Dpsi"]), diff(psi)/dx, tol=tol[2]),# psi'  == psip
      all.equal(m.psi[xn0,"wgt"], (psi/x)[xn0], tol= tol[1]/10))# psi/x == wgt
}
