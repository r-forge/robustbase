##-*- mode: R; kept-new-versions: 50; kept-old-versions: 50 -*-

#### MT Estimators:  [M]-Estimators based on [T]ransformations
#### -------------   Valdora & Yohai (2013)


##' Defining the spline to compute the center of the rho function
##' @title Provide  mu(lambda) as spline function
##' @param cw tuning parameter for rho
##' @return a function, resulting from \code{\link{splinefun}}
##' @author Martin Maechler
mk.m_rho <- function(cw, sFile = paste0("MTesSpl_", format(cw), ".rda")) {
## FIXME: Solution without files, but rather cache inside an environment
## ------  For the default  cw, cache even in robustbase namespace!
## Instead of saving splinefun() ... just save  (la, mm.la), it is much smaller
    useFile <- file.exists(sFile)
    if (useFile) { ## load the spline
        load(sFile)#-> 'cw.0' and 'maprox'
        cw.ok <- (abs(cw - cw.0) < 0.001) # if the last cw was very close..
    }
    if(!useFile || !cw.ok) {
        la <- c(seq(0,2.9,0.1), seq(3,100))# <- MM: I plan to improve this {dependent on cw!}
        mm.la <- rep(0,length(la))

        for(i in 1:length(la))
        { mm.la[i] <- optim(sqrt(la[i]),espRho,method = "L-BFGS-B",lam = la[i],cw = cw)$par }

        maprox <- splinefun(la, mm.la, method = "monoH.FC")
        cw.0 <- cw
        save(cw.0, maprox, file = sFile)
    }
    uu <- maprox
    uu
}
## will be passed around as 'uu'  --> and is used in  mm(.)  below

#######################################################





rho <- function(x,cw) {
    ## defining the rho function
    ans <- rep(0, length(x))
    for (i in 1:length(x)) {
        ans[i] <- min(1-(1-(x[i]/cw)^2)^3,1)
    }
    ans
}


##############################################################

espRho <- function(lam, xx, cw)
{
## compute  E_lambda [ rho_{cw}( sqrt(Y)-xx ) ],  given  (lambda, xx, cw)
    suby <- seq(trunc((max(0,xx-cw))^2),
                trunc((xx+cw)^2)+1)
    subsuby <- suby[rho(sqrt(suby)-xx,cw) < 1]
    terminos <- rho(sqrt(subsuby)-xx,cw)*dpois(subsuby,lam)
    long <- length(subsuby)
    primero <- subsuby[1]
    ultimo <- subsuby[long]
    if (long > 0) ans <- sum(terminos)+1-ppois(ultimo,lam)+ppois(primero-1,lam) else ans <- 1
    ans
}

#################################################

##' @title Compute  m(lambda) := the value of x minimizing espRho(lambda, x, cw)

##' @param lam numeric vector of non-negative values \lambda
##' @param uu the spline function to be used for lambda <= 100, from mk.m_rho()
##' @return
mm <- function(lam,uu)
{
    z <- (lam < 100)
    m <- z*uu(lam)+(1-z)*sqrt(lam)
    m
}## MM: FIXME  mm(Inf, *) currently gives NaN, as 0*Inf = NaN (instead of Inf !)


###############################################################################


weightsMT <- function(x,iw)
{
    p <- ncol(x)
    n <- nrow(x)
    w <- rep(1,n)
    if(iw == 1)
    {
        w <- rep(0,n)
        require(rrcov)
        sest <- CovSest(x[,2:p],method = "bisquare")
        mu <- getCenter(sest)
        sigma <- getCov(sest)
        insigma <- solve(sigma)

        for ( i in 1:n)
        {
            z <- as.vector(x[i,2:p])

            w[i] <- ((z-mu)%*%insigma%*%(z-mu) <= qchisq(0.975,p-1))
        }
    }
    w
}

#######################################################################

##' @title Compute the loss function for MT-glmrob()
##' @param bt beta  (p - vector)
##' @param x  design matrix (n x p)
##' @param y  (Poisson) response  (n - vector)
##' @param cw tuning parameter 'c' for rho_c(.)
##' @param w weight vector (length n)
##' @param uu the spline for the inner part of  m_c(.)
##' @return \sum_{i=1}^n  w_i  \rho(\sqrt(y_i) - m( g(x_i ` \beta) ) )
##'       where  g(.) = exp(.) for the Poisson family
sumaConPesos <- function(bt,x,y,cw, w,uu) {
    n <- nrow(x)
    s <- rep(0,n)
    for (i in 1:n) {
        if (abs(t(bt)%*%x[i,]) < 500)
            s[i] <- rho(sqrt(y[i])-mm(exp(t(bt)%*%x[i,]),uu), cw) else s[i] <- 1
    }
    sum(s*w)
}

###############################################################################

beta0IniCP <- function(x,y,cw,w, uu, nsubm, trace.lev = 1)
{
    ## computes  the initial estmate usng subsamplng with concentration step
    p <- ncol(x)
    n <- nrow(x)
    dev <- matrix(rep(0,n))
    estim.ini <- matrix(0,nsubm,p)
    s2 <- rep(0,nsubm)
    kk <- 0

    for(l in 1:nsubm)
    {
        if(trace.lev) cat(".", if(l %% 50 == 0)paste0(" ",l,"\n"))

        submaux <- sample(n,p )
        estim0 <- betaExacto(submaux,x,y)
        estim0 <- as.vector(estim0)
        correc <- y == 0
        ymueAux <- y+correc
        ## TODO eta <- as.vector(x %*% estim0)
        dev <- abs(y*(log(ymueAux)-c(estim0%*%t(x)))-(y-exp(c(estim0%*%t(x)))))
        devOr <- sort(dev)# FIXME: use  sort(partial = half)
        podador <- dev <= devOr[ceiling(n/2)] # those smaller-equal than lo-median(.)
        xPod <- x[podador,]
        yPod <- y[podador]
	fitE <- tryCatch(glm(yPod ~ xPod-1, family = poisson()),
			 error = function(e)e)
        if(inherits(fitE, "error")) {
            message("glm(.) {inner subsample} error: ", fitE$message)
            s2[l] <- Inf
        } else { ## glm() call succeeded
            betapod <- fitE$coefficients
            estim.ini[l,] <- betapod
### FIXME: replacing  NA's by 0
            ## 1. can be made much more efficiently (notably for *no* NAs)
            ## 2. MM: why is this a good idea???  Rather drop it and try again ?!!
### FIXME this is ~=  betaExacto() .. rather use that here as well
            sinNas <- na.exclude(betapod)
            long <- length(sinNas)
            if(long == p)
                kk <- kk+1
            lugaresNas <- na.action(sinNas)[1:(p-long)]
            estim.ini[l,lugaresNas] <- 0
            ## FIXME: rather only keep the best sample !! --save storage

            s2[l] <- sumaConPesos(estim.ini[l,],x,y,cw, w,uu)
        }
    }

    s0 <- order(s2)
    beta0ini <- estim.ini[s0[1],]

    list(beta = beta0ini, nOksamples = kk)
}## beta0IniCP()

#####################################################################
betaExacto <- function(sbmn,x,y)
{
    ## to each subsample assign the maximum likelihood estimator and
    ## fixing the case mle has NA components
    nmue <- nrow(x)
    p <- ncol(x)
    fitE <- tryCatch(glm(y[sbmn]~x[sbmn,]-1,family = poisson()),
                     error = function(e)e)
    if(inherits(fitE, "error")) {
        message("betaExacto glm(.) error: ", fitE$message)
        return(rep(NA_real_, p))
    }
    ## else --  glm() succeeded
    betaexacto <- fitE $ coefficients
    sinNas <- na.exclude(betaexacto)
    long <- length(sinNas)
    lugaresNas <- na.action(sinNas)[1:(p-long)]
    betaexactoSinNas <- betaexacto
    betaexactoSinNas[lugaresNas] <- 0
    betaexactoSinNas
}



###################################################################################
glmrobMT <- function(x,y, cw = 2.1, iw = 0, nsubm = 500, maxitOpt = 200, tolOpt = 1e-6, trace.lev = 1)
{
    ## MAINFUNCTION  Computes the MT or WMT estimator  for Poisson regression  with intercept starting from the estimator computed in the function
    ## beta0IniC.
    ## INPUT
    ## x  design   matrix with nrows and p columns.
    ## y respone  vector of lengtg n
    ## cw tuning constant. Default value 2.1
    ## iweigths indicator for weights penalizing high leverage points, iweights=1 indicates to use weights iweights=0
    ## indicate notto use way. Default value is iw=0, Our simulation study suggests not to use weights.
    ## nsubm Number of subsamples. Default calue nsubm=500
    ## OUTPUT
    ##$initial is the inital estimate (first component is the intercept)
    ##$final is the final estimate (first component is the intercept)
    ##$nsamples is the number of well  conditioned  subsamples
    ## REQUIRED PACKAGES: tools, rrcov
    uu <- mk.m_rho(cw)
    n <- nrow(x)
    x1 <- cbind(rep(1,n),x)
    w <- weightsMT(x1,iw)
    if(trace.lev) cat("Computing initial estimate with ", nsubm, " sub samples:\n")
    out <- beta0IniCP(x1, y, cw = cw, w = w, uu = uu, nsubm = nsubm, trace.lev = trace.lev)
    beta0 <- out[[1]]

    oCtrl <- list(trace = trace.lev, maxit = maxitOpt,
                  ## "L-BFGS-B" specific
                  lmm = 9, factr = 1/(10*tolOpt))

    if(trace.lev) cat("Optim()izing  sumaConPesos()\n")
### FIXME: see very slow convergence e.g. for the Possum data example
### -----  maybe improve by providing gradient ??
    estim2 <- optim(beta0, sumaConPesos, method = "L-BFGS-B",
                    x = x1, y = y, cw = cw, w = w, uu = uu, control = oCtrl)
    if(estim2$convergence) ## there was a problem
        warning("optim(.) non-convergence: ", estim2$convergence,
                if(nzchar(estim2$message)) paste0("\n", estim2$message))

    out <- list("initial" = beta0, "final" = estim2$par, nsubm = nsubm, "nOksub" = out[[2]],
                converged = (estim2$convergence == 0), optim.counts = estim2$counts,
                optim.control = oCtrl, iw=iw, cw=cw, weights = w)
    out
}
