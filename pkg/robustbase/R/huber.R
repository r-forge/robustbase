
##  A modified "safe" (and more general) Huber estimator:
huberM <-
    function(x, k = 1.5, weights = NULL,
	     tol = 1e-06,
	     mu = if(is.null(weights)) median(x) else wgt.himedian(x, weights),
	     s = if(is.null(weights)) mad(x, center=mu)
		 else wgt.himedian(abs(x - mu), weights),
	     se = FALSE,
	     warn0scale = getOption("verbose"))
{
    ## Author: Martin Maechler, Date: 6 Jan 2003, ff

    ## implicit 'na.rm = TRUE':
    if(any(i <- is.na(x))) {
        x <- x[!i]
        if(!is.null(weights)) weights <- weights[!i]
    }
    n <- length(x)
    sum.w <-
        if(!is.null(weights)) {
            stopifnot(is.numeric(weights), weights >= 0, length(weights) == n)
            sum(weights)
        } else n
    it <- 0L
    NA. <- NA_real_
    if(sum.w == 0) # e.g 'x' was all NA
	return(list(mu = NA., s = NA., it = it, se = NA.)) # instead of error

    if(se && !is.null(weights))
	stop("Std.error computation not yet available for the case of 'weights'")
    missS <- missing(s)
    if(missS && is.na(s)) # e.g. for x = c(-Inf, 1)
        s <- 0
    if (s <= 0) {
        if(s < 0) stop("negative scale 's'")
        if(warn0scale && n > 1)
            warning("scale 's' is zero -- returning initial 'mu'")
    }
    else {
        wsum <- if(is.null(weights)) sum else function(u) sum(u * weights)
	repeat {
	    it <- it + 1L
            y <- pmin(pmax(mu - k * s, x), mu + k * s)
	    mu1 <- wsum(y) / sum.w
	    if (is.na(mu1) || abs(mu - mu1) < tol * s)
		break
	    mu <- mu1
	}
    }
    list(mu = mu, s = s, it = it,
         SE = if(se) s * sqrt(tauHuber(x, mu=mu, s=s, k=k) / n) else NA.)
}

## this is a compatible improvement of MASS' huber() :
## 1) returning median() if mad()=0
## 2)	"	NA when y has only NAs (or length 0)

if(FALSE)
huber <- function (y, k = 1.5, tol = 1e-06)
{
    y <- y[!is.na(y)]
    n <- length(y)
    if(n == 0) # e.g 'y' was all na
	return(list(mu = NA, s = NA))# instead of error
    mu <- median(y)
    s <- mad(y, center=mu)
    if (s == 0) { # FIXME?  make this warning optional
	if(n > 1) warning("scale MAD is zero for this sample")
    }
    else repeat {
	yy <- pmin(pmax(mu - k * s, y), mu + k * s)
	mu1 <- sum(yy)/n
	if (abs(mu - mu1) < tol * s)
	    break
	mu <- mu1
    }
    list(mu = mu, s = s)
}


## Originally from  /u/ftp/NDK/Source-NDK-9/R/rg2-fkt.R :
tauHuber <- function(x, mu, k=1.5, s = mad(x), resid = (x - mu)/s) {
  ## Purpose: Correction factor Tau for the variance of Huber-M-Estimators
  ## -------------------------------------------------------------------------
  ## Arguments: x = data, mu = location, k = tuning parameter of Huber Psi-function
  ## -------------------------------------------------------------------------
  ## Author: Rene Locher Update: R. Frisullo 23.4.02;  M.Maechler (as.log(); s, resid)
  inr <- abs(resid) <= k
  ## psi  <- ifelse(inr, resid, sign(resid)*k)             # psi (x) -- more efficiently:
  psi <- resid; out <- which(!inr); if(length(out)) psi[out] <- sign(resid[out]) * k
  psiP <- as.logical(inr)# = ifelse(abs(resid) <= k, 1, 0) # psi'(x)
  length(x) * sum(psi^2) / sum(psiP)^2
}

