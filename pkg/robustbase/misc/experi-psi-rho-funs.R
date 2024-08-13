#### These are experiments, using the definitions in
#### ==================
#### ../R/psi-rho-funs.R  <<<<<<<< see there
#### ===================


## NB:  "Recycle" code in
## -- /u/maechler/R/MM/STATISTICS/robust/psi-funs.R
##      ~~~~~~~~~~  and
##    /u/maechler/R/MM/STATISTICS/robust/wt-funs.R
##      ~~~~~~~~~~

Hf0 <- function(x, c=1.35) pmax.int(-c, pmin.int(x,c))
Hf <- new("functionX", Hf0)
stopifnot(validObject(Hf)) # ok !


### psiFunc() examples
## classical {trivial, not interesting}:
F1 <- function(x) rep.int(1, length(x))
cPsi <- psiFunc(rho = function(x) x^2 / 2, psi = function(x) x,
                wgt = F1, Dpsi = F1,
                Erho = function(x) rep.int(1/2, length(x)),
                Epsi2 = F1, EDpsi = F1)


### MASS  -- ?rlm --- has
##
##-      psi.huber   (u, k = 1.345, deriv = 0)
##-      psi.hampel  (u, a = 2, b = 4, c = 8, deriv = 0)
##-      psi.bisquare(u, c = 4.685, deriv = 0)
##   where deriv = 0 :  psi(x)/x i.e.  'wgt'

## MM has more in psi-funs.R (see above)



## Reproduce Table 1, p.138 of Hampel, Rousseeuw, Ronchetti, Stahel (1986):
##           ---------------
b <- c(seq(0, 3, by = 0.1), 4, 5, Inf)
A <- huberPsi@Epsi2(b)
B <- huberPsi@EDpsi(b)

options(width = 100, digits = 7)
cbind(b = b, A = A, B = B, V = A/B^2, e = B^2/A,
      gamma. = b/B, k. = 1 + b^2/A, lambda. = 1/B)


(huberPsi2 <- chgDefaults(huberPsi, k = 2)) # works too!

## TODO:  Hampel,  Biweight
## Now the default in robustbase -- .Mpsi.tuning.default("hampel")
Ha.. <- chgDefaults(hampelPsi, k = c(1.5, 3.5, 8) * 0.9016085)
(asVar.. <- Ha..@Epsi2() / Ha..@EDpsi()^2) # 1.052602
1/asVar.. #  0.950027 .. well ..- ==> find a better factor than 0.9016085:

effHa <- function(f) {
    stopifnot(is.finite(f), length(f) == 1L, f > 0)
    cc <- f * c(1.5, 3.5, 8)
    ## return asymp. relative efficiency = 1 / asymp.Var = B^2 / A
    ## eff =
    hampelPsi@EDpsi(cc)^2 / hampelPsi@Epsi2(cc)
}

effHa(0.9016085) #  0.950027
stopifnot(all.equal(0.95, effHa(0.9016085), tolerance = 5e-5)) # 2.8388e-5

findFHa <- function(eff = 0.95, ...) {
    stopifnot(is.finite(eff), length(eff) == 1L, 0 < eff, eff <= 1)
    ur <- uniroot(function(f) effHa(f) - eff, ...)
}
str(r1 <- findFHa(0.95, interval = c(0.88, 0.92)), digits.d = 8)
 ## $ root      : num 0.90144316
 ## $ f.root    : num -1.0108246e-07
 ## $ iter      : int 3
 ## $ init.it   : int NA
 ## $ estim.prec: num 6.1035156e-05
str(r2 <- findFHa(0.95, interval = c(0.88, 0.92), tol = 1e-12), digits.d = 12)
## $ root      : num 0.901443781864
## $ f.root    : num 7.77156117238e-16
## $ iter      : int 5
## $ init.it   : int NA
## $ estim.prec: num 5.00377517199e-13
str(rX <- findFHa(0.95, interval = c(0.88, 0.92), tol = 1e-100, maxiter = 11), digits.d = 16)
##                   vvvvvvvvvvvvvvvvvv  the "correct" factor !!
## $ root      : num 0.9014437818636579
## $ f.root    : num 1.11..e-16
## $ iter      : int 8
## $ init.it   : int NA
## $ estim.prec: num 4.44..e-16

Ha..Corr <- chgDefaults(hampelPsi, k = c(1.5, 3.5, 8) *  0.9014437818636579)
(asVar.. <- Ha..Corr@Epsi2() / Ha..Corr@EDpsi()^2) # 1.052632
1/asVar.. # 0.95
## rel.error:
(1/asVar..)/0.95 -1 # 2.22e-16 (1-2 bits only!)

## Same for Huber:
effHub <- function(k) {
    stopifnot(is.finite(k), length(k) == 1L, k > 0)
    ## return asymp. relative efficiency = 1 / asymp.Var = B^2 / A
    ## eff =
    huberPsi@EDpsi(k)^2 / huberPsi@Epsi2(k)
}

effHub(k = 1.345) # k = 1.345 "literature" default for eff = 0.95
## 0.9500003  ... is indeed, very accurate
stopifnot(all.equal(0.95, effHub(1.345), tolerance = 5e-7)) # need 2.733714e-07

findFHub <- function(eff = 0.95, ...) {
    stopifnot(is.finite(eff), length(eff) == 1L, 0 < eff, eff <= 1)
    ur <- uniroot(function(f) effHub(f) - eff, ...)
}
str(rXhub <- findFHub(0.95, interval = c(1, 1.5), tol = 1e-100, maxiter = 11), digits.d = 16)
## $ root      : num 1.344997508512041
## $ f.root    : num 0
## $ iter      : int 8
## $ init.it   : int NA
## $ estim.prec: num 3.774758e-15



## Experiments only:--- delete later:
PL <- list(rho =
            function(x, k) {
### UNFINISHED
                u <- abs(x)
                lrg <- k[3] <= u
                I <- u < k[1]
                mid <- !I & !lrg    # contains constant and descending
                r <- numeric(length(x))
                r[ I] <- u[I]^2 / 2
                r[mid] <- k*(u[!I] - k / 2)
                r
            },
            psi = function(x, k) {
                ## this is "optimized" ==> factors faster than using ifelse()!
                u <- abs(x)
                lrg <- k[3] <= u
                mid <- k[1] < u & !lrg # constant _and_ descending
                ## x is result for |x| < k[1]
                x[lrg] <- 0
                if(any(mid))
                    x[mid] <- k[1] * sign(x[mid])*
                        pmin.int(1, (u[mid] - k[3])/(k[2] - k[3]))
                x
            },
            wgt  = function(x, k) {
                x <- abs(x)
                lrg <- k[3] <= x
                I <- x < k[1]
                mid <- !I & !lrg    # contains constant and descending
                x[I] <- 1
                x[lrg] <- 0
                if(any(mid))
                    x[mid] <- k[1] / x[mid] *
                        pmin.int(1, (x[mid] - k[3])/(k[2] - k[3]))
                x
            }
           )
plot(function(x) PL$psi(x, k = c(2,4,8)), -10, 10)
## Compare psi(x) / x  with wgt(x) :
plot(function(x) PL$psi(x, k = c(2,4,8))/x, -10, 10)
curve(PL$wgt(x, k = c(2,4,8)), add = TRUE, col="pink")
abline(v = outer(c(-1,1),c(2,4,8)), col="light blue", lty = "1F")
## -> seems fine

ppsi <- PL$psi; formals(ppsi)[-1] <- list(k = c(2,4,8))
str(ppsi)
plot(ppsi, -10, 10) # all looks fine ...
pp <- new("functionX", ppsi)
## now ok

pwgt <- PL$wgt; formals(pwgt)[-1] <- list(k = c(2,4,8))
str(pwgt)
plot(pwgt, -10, 10) # all looks fine ...
pp <- new("functionX", pwgt)## now ok

expre
prho <- PL$rho; formals(prho)[-1] <- list(k = c(2,4,8))
str(prho)
plot(prho, -10, 10) # all looks fine ...
pp <- new("functionX", prho)## now ok

###--- Compute  E[rho]  _numerically_ for Huber:
rho <- function (x, k)
{
    r <- u <- abs(x)
    I <- u < k
    r[I] <- u[I]^2/2
    r[!I] <- k * (u[!I] - k/2)
    r
}

kk <- c(seq(0.01, 0.5 - 1/64, by=1/64), seq(0.5, 1.9875, by = 1/32), seq(2, 4, by = 1/16))
ErhoA <- lapply(kk, function(k)
                integrate(function(x)rho(x, k=k) * dnorm(x), -Inf,Inf,
                          rel.tol = 1e-14))
str(Erho <- sapply(ErhoA, "[[", "value"))
plot(Erho ~ kk, type = 'o', col = 2, cex = 1/2)
## Compare with "exact":
lines(kk, huberPsi@Erho(kk), col = adjustcolor("forestgreen", 1/3), lwd=3)
## seems very accurate  ... but then we can only see errors " O(1/1000) "

plot(kk, relErrV(huberPsi@Erho(kk), Erho), type="l",
     main = quote("relative Error of integrate"*{}(rho[list(Huber, k)])), xlab = quote(k))
abline(h = 0, lty=3)
## ok: errors ~=  rel.tol = 1e-14  from integrate() above



## Hampel: in {robustbase} now, ---> ../R/psi-rho-funs.R
## ------
hampelPsi # {short print() method !}

## TODO:  Biweight :
if(FALSE)
tukeyPsi <- c() ##########

## to use for other types, just change 'ggw' argument
## standardized to have Dpsi(0) = 1
## to have rho(inf) = 1 use .M.chi instead (as well as deriv + 1)
## using this results in an error while preparing for lazy loading:
## (MM, MK: the error arises during the validity check)
## ** preparing package for lazy loading
## Creating a generic function from function "chgDefaults"
## Error in .M.psi(x, k, "ggw", -1) : object 'R_psifun' not found
## Error : unable to load R code in package 'robustbase'
## ERROR: lazy loading failed for package ‘robustbase’
## ('R_psifun' is the pointer to the C-function used in .M.psi)
ggwPsi <- psiFunc(rho = function(x, k) .M.psi(x, k, 'ggw', -1),
                  psi = function(x, k) .M.psi(x, k, 'ggw', 0),
                  wgt = function(x, k) .M.wgt(x, k, 'ggw'),
                  Dpsi = function(x, k) .M.psi(x, k, 'ggw', 1),
                  Erho = function(k) lmrob.E(.M.psi(r, k, 'ggw', -1),
                    use.integrate = TRUE),
                  Epsi2 = function(k) lmrob.E(.M.psi(r, k, 'ggw', 0)^2,
                    use.integrate = TRUE),
                  EDpsi = function(k) lmrob.E(.M.psi(r, k, 'ggw', 1),
                    use.integrate = TRUE),
                  k = c(-0.5, 1.5, 0.95, NA))

## maybe TODO: Optimal tanh() estimator for location


### --> scale - rho/psi functions
