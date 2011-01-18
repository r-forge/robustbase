require(robustbase)

## Demonstrate that  one of  tukeyChi() / tukeyPsi1() is superfluous

EQ <- function(x,y) all.equal(x,y, tol = 1e-13)

x <- seq(-4,4, length=201)
c. <- pi

stopifnot(EQ(tukeyChi(x, c.),
             6/c.^2* tukeyPsi1(x, c., deriv=-1)),
          EQ(tukeyChi(x, c., deriv= 1),
             6/c.^2* tukeyPsi1(x, c., deriv= 0)),
          EQ(tukeyChi(x, c., deriv= 2),
             6/c.^2* tukeyPsi1(x, c., deriv= 1)))

## Test if default arguments are used
tPsi <- chgDefaults(huberPsi, k = 2)

stopifnot(tPsi@rho(1:10, k=2) == tPsi@rho(1:10),
          tPsi@psi(1:10, k=2) == tPsi@psi(1:10),
          tPsi@Dpsi(1:10, k=2) == tPsi@Dpsi(1:10),
          tPsi@wgt(1:10, k=2) == tPsi@wgt(1:10))

## FIXME: Default arguments are not used for E... slots
try(tPsi@Erho())
try(tPsi@Epsi2())
try(tPsi@EDpsi())
