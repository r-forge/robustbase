##-*- mode: R; kept-new-versions: 50; kept-old-versions: 50 -*-

#### MT Estimators:  [M]-Estimators based on [T]ransformations
#### -------------   Valdora & Yohai (2013)


##' Defining the spline to compute the center of the rho function
##' @title Provide  mu(lambda) as spline function
##' @param cw tuning parameter for rho
##' @return a function, resulting from \code{\link{splinefun}}
##' @author Victor Yohai;  many changes: Martin Maechler
mk.m_rho <- function(cw,
 opt.method = c("L-BFGS-B", "Brent", "Nelder-Mead", "BFGS", "CG", "SANN"),
##optim(): method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
## MM: 'Brent' seems best overall
                     lambda = c(seq(0,2.9, by=0.1), seq(3,100)),
                     reltol = sqrt(.Machine$double.eps), trace = 0,
		     sFile = paste0("MTesSpl_", format(cw), ".rda"),
                     recompute  = getOption("robustbase:m_rho_recompute", FALSE))
{
## FIXME: Solution without files, but rather cache inside an environment
## ------  For the default  cw, cache even in robustbase namespace!
## Instead of saving splinefun() ... just save	(lambda, mm.la), it is much smaller

    if(recompute) {
	useFile <- FALSE
    } else {
	useFile <- file.exists(sFile)
	if (useFile) { ## load the spline
	    load(sFile)#-> 'cw.0' and 'm.approx'
	    ## check if its cw was very close to this one:
	    if(cw.ok <- is.numeric(cw0 <- environment(m.approx)$cw))
		cw.ok <- (abs(cw - environment(m.approx)$cw) < 0.001)
	}
    }
    if(!useFile || !cw.ok) {
	nl <- length(lambda)
	mm.la <- numeric(nl)
	s.la <- sqrt(lambda)

        ## MM: Speedwise,   Brent > L-BFGS-B > BFGS > ..  for cw >= ~ 1.5
        ##         L-BFGS-B > Brent  for  cw ~= 1
	opt.method <- match.arg(opt.method)
	oCtrl <- list(reltol=reltol, trace=trace)
	if(opt.method %in% c("Brent", "L-BFGS-B")) { ## use bounds
	    if(opt.method == "L-BFGS-B")# yuck! why is this necessary!!
		oCtrl <- list(factr = 1/(10*reltol), trace=trace)
	    for(i in seq_len(nl))
		mm.la[i] <- optim(s.la[i], espRho, lam = lambda[i], cw = cw,
				  method = opt.method, control = oCtrl,
				  lower = 0, upper = .01 + 2*s.la[i])$par
	} else {
	    for(i in seq_len(nl))
		mm.la[i] <- optim(s.la[i], espRho, lam = lambda[i], cw = cw,
				  method = opt.method, control = oCtrl)$par
	}
	m.approx <- splinefun(lambda, mm.la, method = "monoH.FC")
        e <- environment(m.approx)
	assign("lambda.max", max(lambda), envir=e)
	assign("cw", cw, envir=e)
	save(m.approx, file = sFile)
    }
    m.approx
}
## result 'm.approx' will be  used in mm(.), and "everywhere" below

#######################################################





##' Tukey's Bisquare (aka "biweight") rho function: rho~() = rho scaled to have rho(Inf) = 1
rho <- function(x,cw) pmin(1, 1 - (1-(x/cw)^2)^3)
## even faster:
rho <- function(x,cw) lmrob.chifun(x, cc=cw, psi="tukey")
## NB: in sumaConPesos(), mm(.), ... we make use of the fact  that  rho(Inf) = 1

##############################################################

espRho <- function(lam, xx, cw)
{
## compute  E_lambda [ rho_{cw}( sqrt(Y)-xx ) ],  given  (lambda, xx, cw)
## for Y ~ Pois(lambda) ;  rho(.) = Tukey's Bisquare
    k <- seq(as.integer((max(0,xx-cw))^2), as.integer((xx+cw)^2)+1L)
    inner <- (rhoS.k <- rho(sqrt(k)-xx, cw)) < 1
    ii <- k[inner]
    terminos <- rhoS.k[inner] * dpois(ii,lam)
    if((len.ii <- length(ii)) > 0) {
        primero <- ii[1]
        ultimo <- ii[len.ii]
        sum(terminos) + 1-ppois(ultimo,lam) + ppois(primero-1,lam)
    } else 1
}

#################################################

##' @title Compute  m(lambda) := the value of x minimizing espRho(lambda, x, cw)

##' @param lam numeric vector of non-negative values \lambda
##' @param m.approx the spline function to be used for "small" lambda, from mk.m_rho()
##' @return
mm <- function(lam, m.approx)
{
    la.max <- environment(m.approx)$lambda.max
    z <- ((m <- lam) <= la.max)
    m[z] <- m.approx(lam[z])
    if(any(i <- !z)) m[i] <- sqrt(lam[i])
    m
}


###############################################################################

##' @title Compute the loss function for MT-glmrob()
##' @param beta beta  (p - vector)
##' @param x  design matrix (n x p)
##' @param y  (Poisson) response  (n - vector)
##' @param cw tuning parameter 'c' for rho_c(.)
##' @param w weight vector (length n)
##' @param m.approx the spline for the inner part of  m_c(.)
##' @return \sum_{i=1}^n  w_i  \rho(\sqrt(y_i) - m( g(x_i ` \beta) ) )
##'       where  g(.) = exp(.) for the Poisson family
sumaConPesos <- function(beta,x,y,w, cw, m.approx) {
    eta <- x %*% beta
    s <- rho(sqrt(y) - mm(exp(eta), m.approx), cw)
    sum(s*w)
}

###############################################################################

beta0IniCP <- function(x,y,cw,w, m.approx, nsubm, trace.lev = 1)
{
    ## computes  the initial estmate usng subsamplng with concentration step
    stopifnot(is.matrix(x), (nsubm <- as.integer(nsubm)) >= 1)
    p <- ncol(x)
    n <- nrow(x)
    s2.best <- Inf; b.best <- rep(NA_real_, p)
    kk <- 0
    for(l in 1:nsubm)
    {
        if(trace.lev) {
            if(trace.lev > 1)
                cat(sprintf("%3d:",l))
            else cat(".", if(l %% 50 == 0) paste0(" ",l,"\n"))
        }
        i.sub <- sample(n,p )
        estim0 <- as.vector( betaExacto(x[i.sub,], y[i.sub]) )
	if(any(is.na(estim0))) ## do not use it
	    next
        eta <- as.vector(x %*% estim0)
        ## adev :=  abs( 1/2 * dev.residuals(.) )
        ## y+(y==0) : log(0) |-> -Inf fine; but if eta == -Inf, we'd get NaN
        adev <- abs(y*(log(y+(y == 0)) - eta) - (y-exp(eta)))
        ## poisson()'s dev.resids():  2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
        ##   == 2*wt* (y * ifelse(y == 0, 0, log(y) - log(mu))  - (y - mu))
        ##   == 2*wt* ifelse(y == 0, mu, y*(log(y) - log(mu)) - (y - mu))
        ##  where mu <- exp(eta)

        if(trace.lev > 1) cat(sprintf(" D=%11.7g ", sum(adev)))
        half <- ceiling(n/2)
        srt.d <- sort(adev, partial=half)
        podador <- adev <= srt.d[half] # those smaller-equal than lo-median(.)
        xPod <- x[podador,]
        yPod <- y[podador]
	fitE <- tryCatch(glm(yPod ~ xPod-1, family = poisson()),
			 error = function(e)e)
        if(inherits(fitE, "error")) {
            message("glm(.) {inner subsample} error: ", fitE$message)
            if(trace.lev > 1) cat("\n")
            ## s2[l] <- Inf
        } else { ## glm() call succeeded
	    betapod <- as.vector( fitE$coefficients )
	    if(any(is.na(betapod))) ## do not use it
		next
	    kk <- kk+1
	    s2 <- sumaConPesos(betapod, x=x, y=y, w=w,
                               cw=cw, m.approx=m.approx)
            ## estim.ini[l,] <- betapod
            if(trace.lev > 1) cat(sprintf("s2=%14.9g", s2))
            if(s2 < s2.best) {
                if(trace.lev > 1) cat(" New best!\n")
                b.best <- betapod
                s2.best <- s2
            } else if(trace.lev > 1) cat("\n")
        }
    }

    ## s0 <- order(s2)
    ## beta0ini <- estim.ini[s0[1],]

    list(beta = b.best, nOksamples = kk, s2 = s2.best)
}## beta0IniCP()

#####################################################################
betaExacto <- function(x,y)
{
    ## to each subsample assign the maximum likelihood estimator and
    ## fixing the case mle has NA components
    p <- ncol(x)
    fitE <- tryCatch(glm(y ~ x-1, family = poisson()),
                     error = function(e)e)
    if(inherits(fitE, "error")) {
        message("betaExacto glm(.) error: ", fitE$message)
        return(rep(NA_real_, p))
    }
    ## else --  glm() succeeded
    if(FALSE) { ## original; MM finds it unneeded
        beta. <- fitE $ coefficients
        sinNas <- na.exclude(beta.)
        long <- length(sinNas)
        lugaresNas <- na.action(sinNas)[1:(p-long)]
        beta.SinNas <- beta.
        beta.SinNas[lugaresNas] <- 0
        beta.SinNas
    }
    fitE $ coefficients
}



###################################################################################
glmrobMT <- function(x,y, cw = 2.1, wts = "none", nsubm = 500,
                     maxitOpt = 200, tolOpt = 1e-6, trace.lev = 1)
{
    ## MAINFUNCTION  Computes the MT or WMT estimator  for Poisson regression  with intercept starting from the estimator computed in the function
    ## beta0IniC.
    ## INPUT
    ## x  design   matrix with nrows and p columns.
    ## y respone  vector of length n
    ## cw tuning constant. Default value 2.1
    ## iweigths indicator for weights penalizing high leverage points, iweights=1 indicates to use weights iweights=0
    ## indicate notto use way. Default value is iw=0, Our simulation study suggests not to use weights.
    ## nsubm Number of subsamples. Default calue nsubm=500
    ## OUTPUT
    ##$initial is the inital estimate (first component is the intercept)
    ##$final is the final estimate (first component is the intercept)
    ##$nsamples is the number of well  conditioned  subsamples
    ## REQUIRED PACKAGES: tools, rrcov
    m.approx <- mk.m_rho(cw)
    n <- nrow(x)
    x1 <- cbind(rep(1,n),x)
    w <- robXweights(wts, x1)
    if(trace.lev) cat("Computing initial estimate with ", nsubm, " sub samples:\n")
    out <- beta0IniCP(x1, y, cw = cw, w = w, m.approx = m.approx, nsubm = nsubm, trace.lev = trace.lev)
    beta0 <- out[[1]]

    oCtrl <- list(trace = trace.lev, maxit = maxitOpt,
                  ## "L-BFGS-B" specific
                  lmm = 9, factr = 1/(10*tolOpt))

    if(trace.lev) cat("Optim()izing  sumaConPesos()\n")
### FIXME: see very slow convergence e.g. for the Possum data example
### -----  maybe improve by providing gradient ??
    estim2 <- optim(beta0, sumaConPesos, method = "L-BFGS-B",
                    x = x1, y = y, w = w, cw = cw, m.approx = m.approx, control = oCtrl)
    if(estim2$convergence) ## there was a problem
        warning("optim(.) non-convergence: ", estim2$convergence,
                if(nzchar(estim2$message)) paste0("\n", estim2$message))

    list("initial" = beta0, "final" = estim2$par, nsubm = nsubm, "nOksub" = out[[2]],
	 converged = (estim2$convergence == 0), optim.counts = estim2$counts,
	 optim.control = oCtrl, wts=wts, cw=cw, weights = w)
}
