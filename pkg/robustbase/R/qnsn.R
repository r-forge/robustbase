### Note: till 2010, a slightly wrong constant = 2.2219 was in use.
### Error reported by Peter Ruckdeschel, U.Bayreuth, 15.Juli 2010
### correct constant == 1 / (sqrt(2) * qnorm(5/8)) == 2.219144
### -- but wrong constant, 2.2219, is already in the the original Fortran qn.f
Qn.corr <- 2.2219 / 2.21914
##
Qn <- function(x, constant = 2.21914, finite.corr = missing(constant))
{
    ## Purpose: Rousseeuw and Croux's  Q_n() robust scale estimator
    ## Author: Martin Maechler, Date: 14 Mar 2002, 10:43
    n <- length(x)
    if(n == 0) return(NA) else if(n == 1) return(0.)

    r <- constant *
        .C(Qn0, as.double(x), n, res = double(1))$res

    if (finite.corr)
	(if (n <= 9) {
	    if      (n == 2)  .399
	    else if (n == 3)  .995
	    else if (n == 4)  .513
	    else if (n == 5)  .845
	    else if (n == 6)  .612
	    else if (n == 7)  .858
	    else if (n == 8)  .670
	    else if (n == 9)  .873
	} else {
	    if (n %% 2) ## n odd
		n / (n + 1.4)
	    else ## n even
		n / (n + 3.8)
        }
         ) * r
    else r
}

if(FALSE) {
## How to correct the 9 finite sample correction factors which were from simulations
## computed with the original wrong constant?
## We need  [const * fc]_{new} == [const * fc]_{old}.
## ==> fc_{new} = fc_old *  (const_{old} / const_{new}) = fc_old * Qn.corr
Qn.fin.9.old <- ## original finite sample corrections for n = 2:9:
    c(.399, .994, .512, .844, .611, .857, .669, .872)
round(Qn.fin.9.old * Qn.corr, 3)
##-> 0.399 0.995 0.513 0.845 0.612 0.858 0.670 0.873
}


Sn <- function(x, constant = 1.1926, finite.corr = missing(constant))
{
    ## Purpose: Rousseeuw and Croux's  S_n() robust scale estimator
    ## Author: Martin Maechler, Date: 14 Mar 2002, 10:43

    n <- length(x)
    if(n == 0) return(NA) else if(n == 1) return(0.)

    r <- constant * .C(Sn0,
                       as.double(x), n,
                       as.integer(!is.unsorted(x)),# is.sorted
                       res = double(1), a2 = double(n))$res
    ## NB: a2[] could be used for confidence intervals and other estimates!
    if (finite.corr) (
	if (n <= 9) {
            c(0.743, # n = 2
              1.851, 0.954,# 3 & 4
              1.351, 0.993,# 5 & 6
              1.198, 1.005,# 7 & 8
              1.131)[n - 1]
	} else if (n %% 2) # n odd, >= 11
            n / (n - 0.9)
        else # n even, >= 10
            1
    ) * r
    else r
}

wgt.himedian <- function(x, weights = rep(1,n))
{
    ## Purpose: weighted hiMedian of x
    ## Author: Martin Maechler, Date: 14 Mar 2002
    n <- length(x <- as.double(x))
    stopifnot(storage.mode(weights) %in% c("integer", "double"))
    if(n != length(weights))
	stop("'weights' must have same length as 'x'")
    ## if(is.integer(weights)) message("using integer weights")
    .C(if(is.integer(weights)) wgt_himed_i else wgt_himed,
       x, n, weights,
       res = double(1))$res
}

## To be used directly as  'scaleFun'  in  'covOGK()' :
s_Qn <- function(x, mu.too = FALSE, ...)
    c(if(mu.too) median(x), Qn(x, ...))

s_Sn <- function(x, mu.too = FALSE, ...)
    c(if(mu.too) median(x), Sn(x, ...))
