#### http://www.econ.kuleuven.be/public/NDBAE06/programs/roblog/ :
####
#### August 06, 2010  2:14 PM 9121 BYlogreg.r.txt == BYlogreg.R (*this* original)
####    May 04, 2005  9:24 AM 6702 BYlogreg.txt   == BYlogreg.R.~2005~
####    May 04, 2005  9:25 AM 6720 WBYlogreg.txt  == WBYlogreg.R


##  Computation of the estimator of Bianco and Yohai (1996) in logistic regression
##  -------------
##  Christophe Croux, Gentiane Haesbroeck
##  (thanks to Kristel Joossens and Valentin Todorov for improving the code) -
##  ==> Now "contains" both the *weighted* and regular, unweighted  BY-estimator
##
##  This program computes the estimator of Bianco and Yohai in
##  logistic regression. By default, an intercept term is included
##  and p parameters are estimated.
##
##  For more details we refer to
##     Croux, C., and Haesbroeck, G. (2003), ``Implementing the Bianco and Yohai estimator for Logistic Regression'',
##     Computational Statistics and Data Analysis, 44, 273-295
##

## Output:
##
## A list with the follwoing components:
## convergence   - TRUE or FFALSE if convergence achieved or not
## objective     - value of the objective function at the minimum
## coef          - estimates for the parameters
## sterror       - standard errors of the parameters (if convergence is TRUE)
##
##Example:
##
## x0 <- matrix(rnorm(100,1))
## y  <- as.numeric(runif(100)>0.5)        #numeric(runif(100)>0.5)
## BYlogreg(x0,y)
##

BYlogreg <- function(x0, y, initwml=TRUE, intercept=TRUE, const=0.5,
                     kmax=1e3, maxhalf=10)
{
    if(!is.numeric(y))
        y <- as.numeric(y)
    if(!is.null(dim(y))) {
        if(ncol(y) != 1)
            stop("y is not onedimensional")
        y <- as.vector(y)
    }
    n <- length(y)

    if(is.data.frame(x0)) {
        x0 <- data.matrix(x0)
    } else if (!is.matrix(x0)) {
        x0 <- matrix(x0, length(x0), 1,
                     dimnames = list(names(x0), deparse(substitute(x0))))
    }
    if(nrow(x0) != n)
        stop("Number of observations in x and y not equal")

    if(intercept) {
        x <- cbind("Intercept" = 1, x0)
    } else { # x0 := x without the  "intercept column"
        x <- x0
        all1 <- apply(x == 1, 2, all)
        if(any(all1))
            x0 <- x[,!all1, drop = FALSE]
        else message("no intercept in the model")
    }

    na.x <- !is.finite(rowSums(x))
    na.y <- !is.finite(y)
    ok <- !(na.x | na.y)
    x <- x[ok, , drop = FALSE]
    y <- y[ok] # y[ok, , drop = FALSE]
    dx <- dim(x)
    n <- dx[1]
    if (n == 0)
        stop("All observations have missing values!")
    p <- ncol(x)

    ## Smallest value of the scale parameter before implosion
    sigmamin <- 1e-4

    ## Computation of the initial value of the optimization process
    gstart <-
        if(initwml) {
            ## hp <- floor(n*(1-0.25))+1
            ##        mcdx <- cov.mcd(x0, quantile.used =hp,method="mcd")
            ##        rdx=sqrt(mahalanobis(x0,center=mcdx$center,cov=mcdx$cov))
            ## mcdx <- CovMcd(x0, alpha=0.75)
            ## rdx  <- sqrt(getDistance(mcdx))
            mcd <- covMcd(x0, alpha=0.75)
            D <- mahalanobis(mcd$X, mcd$center, mcd$cov)
            D <- sqrt(D)
            vc  <- sqrt(qchisq(0.975, p-1))
            wrd <- D <= vc
            glm(y~x0, family=binomial, subset=wrd)$coef
        } else {
            glm(y~x0, family=binomial)$coef
        }

    sigmastart <- 1/sqrt(sum(gstart^2))
    xistart <- gstart*sigmastart
    stscores <- x %*% xistart
    sigma1 <- sigmastart

    ## Initial value for the objective function
    oldobj <- mean(phiBY3(stscores/sigmastart,y,const))
    kstep <- jhalf <- 1

    while(kstep < kmax & jhalf < maxhalf)
    {
        unisig <- function(sigma) mean(phiBY3(stscores/sigma,y,const))

        optimsig <- nlminb(sigma1, unisig, lower=0)# "FIXME" arguments to nlminb()
        sigma1 <- optimsig$par

        if(sigma1 < sigmamin) {
            warning("Explosion")
            kstep <- kmax
        } else {
            gamma1 <- xistart/sigma1
            scores <- stscores/sigma1
            newobj <- mean(phiBY3(scores,y,const))
            oldobj <- newobj
            gradBY3 <- colMeans((derphiBY3(scores,y,const)%*%matrix(1,ncol=p))*x)
            h <- -gradBY3+ (gradBY3 %*% xistart) *xistart
            finalstep <- h/sqrt(sum(h^2))
            xi1 <- xistart+finalstep
            xi1 <- xi1/sum(xi1^2)
            scores1 <- (x%*%xi1)/sigma1
            newobj <- mean(phiBY3(scores1,y,const))

### stephalving
            hstep <- jhalf <- 1
            while(jhalf <= maxhalf & newobj > oldobj)
            {
                hstep <- hstep/2
                xi1 <- xistart+finalstep*hstep
                xi1 <- xi1/sqrt(sum(xi1^2))
                scores1 <- x%*%xi1/sigma1
                newobj <- mean(phiBY3(scores1,y,const))
                jhalf <- jhalf+1
            }

            if(jhalf == maxhalf+1 & newobj > oldobj)
            {
                print("Convergence Achieved")
            } else {
                jhalf <- 1
                xistart <- xi1
                oldobj <- newobj
                stscores <- x%*% xi1
                kstep <- kstep+1
            }
        }
    }

    if(kstep == kmax) {
        warning("No convergence in ", kstep, " steps.")
        list(convergence=FALSE, objective=0, coef= rep(NA,p))
    } else {
        gammaest <- xistart/sigma1
        list(convergence=TRUE, objective=oldobj, coef=gammaest,
             sterror = sterby3(x0, y, const, gammaest),
             iter = kstep)
    }
}


### Functions needed for the computation of estimator of Bianco and Yohai ----------------------

## From their paper:

## A last remark is worth mentioning: when huge outliers occur in
## the logistic regression setting, often numerical imprecision occurs in the computation
## of the deviances given by
##    d(s;y_i)= -y_i log F(s) - (1-y_i) log{1-F(s)} .
##
## Instead of directly computing this expression, it can be seen that a
## numerically more stable and accurate formula is given by
##    log(1 + exp(-abs(s))) + abs(s)* ((y-0.5)*s < 0)
## in which the second term equals abs(s) if the observation is misclassified, 0 otherwise.
dev1 <- function(s,y) log(1+exp(-abs(s))) + abs(s)*((y-0.5)*s<0)
dev2 <- function(s,y) log1p(exp(-abs(s))) + abs(s)*((y-0.5)*s<0)
dev3 <- function(s,y) -( y  * plogis(s, log.p=TRUE) +
                        (1-y)*plogis(s, lower.tail=FALSE, log.p=TRUE))

log1pexp <- function(x) {
    if(has.na <- any(ina <- is.na(x))) {
	y <- x
	x <- x[ok <- !ina]
    }
    t1 <- x <= 18
    t2 <- !t1 & (tt <- x <= 33.3)
    r <- x
    r[ t1] <- log1p(exp(x[t1]))
    r[ t2] <- { x2 <- x[t2]; x2 + exp(-x2) }
    r[!tt] <- x[!tt]
    if(has.na) { y[ok] <- r ; y } else r
}


phiBY3 <- function(s,y,c3)
{
  s <- as.double(s)
  ## MM FIXME  log(1 + exp(-.)) ... but read the note above !! ---
  dev <- log(1+exp(-abs(s))) + abs(s)*((y-0.5)*s<0)
  rhoBY3(dev,c3) + GBY3Fs(s,c3) + GBY3Fsm(s,c3)
}

rhoBY3 <- function(t,c3)
{
  (t*exp(-sqrt(c3))*as.numeric(t <= c3))+
    (((exp(-sqrt(c3))*(2+(2*sqrt(c3))+c3))-(2*exp(-sqrt(t))*(1+sqrt(t))))*as.numeric(t >c3))
}

psiBY3 <- function(t,c3)
{
    (exp(-sqrt(c3))*as.numeric(t <= c3))+(exp(-sqrt(t))*as.numeric(t >c3))
}

derpsiBY3 <- function(t,c3)
{
    ## MM FIXME: pre-alloc
    res <- NULL
    for(i in seq_along(t)) {
        if(t[i] <= c3)
        {
            res <- rbind(res,0)
        } else
        {
            res <- rbind(res,-exp(-sqrt(t[i]))/(2*sqrt(t[i])))
        }
    }
    res
}

sigmaBY3 <- function(sigma,s,y,c3)
{
    mean(phiBY3(s/sigma,y,c3))
}

derphiBY3 <- function(s,y,c3)
{
    Fs <- exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
    ds <- Fs*(1-Fs)
    dev <- log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
    Gprim1 <- log(1+exp(-abs(s)))+abs(s)*(s<0)
    Gprim2 <- log(1+exp(-abs(s)))+abs(s)*(s>0)

    return(-psiBY3(dev,c3)*(y-Fs)+((psiBY3(Gprim1,c3)-psiBY3(Gprim2,c3))*ds))
}

der2phiBY3 <- function(s, y, c3)
{
    s <- as.double(s)
    Fs <- exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
    ds <- Fs*(1-Fs)
    dev <- log(1+exp(-abs(s)))+abs(s)*((y-0.5)*s<0)
    Gprim1 <- log(1+exp(-abs(s)))+abs(s)*(s<0)
    Gprim2 <- log(1+exp(-abs(s)))+abs(s)*(s>0)
    der2 <- (derpsiBY3(dev,c3)*(Fs-y)^2)+(ds*psiBY3(dev,c3))
    der2 <- der2+(ds*(1-2*Fs)*(psiBY3(Gprim1,c3)-psiBY3(Gprim2,c3)))

    der2 - ds*(derpsiBY3(Gprim1,c3)*(1-Fs) +
               derpsiBY3(Gprim2,c3)*  Fs )
}


GBY3Fs <- function(s,c3)
{
    Fs <- exp(-(log(1+exp(-abs(s)))+abs(s)*(s<0)))
    resGinf <- exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fs))))-1)
    resGinf <- (resGinf+(Fs*exp(-sqrt(-log(Fs)))))*as.numeric(s <= -log(exp(c3)-1))
    resGsup <- ((Fs*exp(-sqrt(c3)))+(exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1)))*as.numeric(s > -log(exp(c3)-1))

    return(resGinf+resGsup)
}


GBY3Fsm <- function(s,c3)
{
    Fsm <- exp(-(log(1+exp(-abs(s)))+abs(s)*(s>0)))
    resGinf <- exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fsm))))-1)
    resGinf <- (resGinf+(Fsm*exp(-sqrt(-log(Fsm)))))*as.numeric(s >= log(exp(c3)-1))
    resGsup <- ((Fsm*exp(-sqrt(c3)))+(exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1)))*as.numeric(s < log(exp(c3)-1))
    resGinf + resGsup
}

## Compute the standard erros of the estimates -
##  this is done by estimating the asymptotic variance of the normal
##  limiting distribution of the BY estimator - as derived in Bianco
##  and Yohai (1996)
##
sterby3 <- function(x0, y, const, estim, intercept=TRUE)
{
    stopifnot(length(dim(x0)) == 2)
    d <- dim(z <- if(intercept) cbind(1, x0) else x0)
    argum <- z %*% estim
    n <- d[1]
    p <- d[2]
    matM <- IFsquar <- matrix(0, nrow=p, ncol=p)
    for(i in 1:n)
    {
        myscalar <- as.numeric(der2phiBY3(argum[i],y[i],const))
        zzt <- tcrossprod(z[i,])
        matM <- matM + myscalar * zzt
        IFsquar <- IFsquar + derphiBY3(argum[i],y[i],const)^2 * zzt
    }

    matM    <- matM/n
    matMinv <- solve(matM)
    IFsquar <- IFsquar/n
    asvBY   <- matMinv %*% IFsquar %*% t(matMinv)

    sqrt(diag(asvBY))/sqrt(n)
}

