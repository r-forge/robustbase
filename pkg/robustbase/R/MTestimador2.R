

spcreation <- function(cw) {
    ##defining the spline to compute the centerof the rho function
    u  <- FALSE
    lf  <- list_files_with_exts(".","mtes",all.files = TRUE,full.names = FALSE)
    if (length(lf) == 2) {
        load("sp.mtes")
        cwguar  <- scan("cw.mtes")
        u  <- (abs(cw-cwguar) < 0.001)
    }

    if(length(lf) != 2 | u == FALSE) {
        la <- c(seq(0,2.9,0.1),seq(3,100))
        mm.la <- rep(0,length(la))

        for(i in 1:length(la))
        { mm.la[i] <- optim(sqrt(la[i]),espRho,method = "L-BFGS-B",lam = la[i],cw = cw)$par }

        maprox <- splinefun(la,mm.la , method = "monoH.FC")
        save("maprox",file = "sp.mtes")
        write(cw,file = "cw.mtes")
    }
    uu  <- maprox
    uu
}

#######################################################





rho <- function(x,cw) {
    ##defining the rho function
    ans <- rep(0, length(x))
    for (i in 1:length(x))
    {
        ans[i] <- min(1-(1-(x[i]/cw)^2)^3,1)
    }
    ans
}


##############################################################

espRho <- function(lam,xx,cw)
{
    ##defining the   function that given  (lambda, xx,cw) computes  E_lambda (rho_{cw}(sqrt(y)-xx))
    suby <- seq(trunc((max(0,xx-cw))^2),trunc((xx+cw)^2)+1)
    subsuby <- suby[rho(sqrt(suby)-xx,cw) < 1]
    terminos <- rho(sqrt(subsuby)-xx,cw)*dpois(subsuby,lam)
    long <- length(subsuby)
    primero <- subsuby[1]
    ultimo <- subsuby[long]
    if (long > 0) ans <- sum(terminos)+1-ppois(ultimo,lam)+ppois(primero-1,lam) else ans <- 1
    ans
}

#################################################


mm <- function(lam,uu)
{
    ##computing  the value of xx minimizing espRho(lambda,xx,cw)
    z  <- (lam < 100)
    m  <- z*uu(lam)+(1-z)*sqrt(lam)
    m
}


###############################################################################


weights <- function(x,iw)
{
    p <- ncol(x)
    n  <- nrow(x)
    w <- rep(1,n)
    if(iw == 1)
    {
        w <- rep(0,n)
        sest <- CovSest(x[,2:p],method = "bisquare")
        mu <- getCenter(sest)
        sigma <- getCov(sest)
        insigma <- solve(sigma)

        for ( i in 1:n)
        {
            z  <- as.vector(x[i,2:p])

            w[i] <- ((z-mu)%*%insigma%*%(z-mu) <= qchisq(0.975,p-1))
        }
    }
    w
}

#######################################################################
sumaConPesos <- function(bt,x,y,cw, w,uu) {
    ##defining  the  function that compute the loss function
    n  <- nrow(x)
    s <- rep(0,n)
    for (j in 1:n)	{

        if (abs(t(bt)%*%x[j,]) < 500) s[j] <- rho(sqrt(y[j])-mm(exp(t(bt)%*%x[j,]),uu),cw) else s[j] <- 1
    }
    sum(s*w)
}

###############################################################################



beta0IniCP <- function(x,y,cw,w,uu,nsubm)
{
    ##defining   the  funcion that computes  the initial estmate usng subsamplng with concentration step
    p <- ncol(x)
    n <- nrow(x)
    dev <- matrix(rep(0,n))
    estim.ini <- matrix(0,nsubm,p)
    s2 <- rep(0,nsubm)
    kk  <- 0

    for(l in 1:nsubm)
    {
        ##print(l)

        submaux <- sample(n,p )
        estim0 <- betaExacto(submaux,x,y)
        estim0 <- as.vector(estim0)
        correc <- y == 0
        ymueAux <- y+correc
        dev <- abs(y*(log(ymueAux)-c(estim0%*%t(x)))-(y-exp(c(estim0%*%t(x)))))
        devOr <- sort(dev)
        podador <- dev <= devOr[ceiling(n/2)]
        xPod <- x[podador,]
        yPod <- y[podador]
        betapod <- glm(yPod ~ xPod-1, family = poisson())$coefficients
        sinNas <- na.exclude(betapod)
        long <- length(sinNas)
        if(long == p)
            kk  <- kk+1
        lugaresNas <- na.action(sinNas)[1:(p-long)]
        estim.ini[l,] <- betapod
        estim.ini[l,lugaresNas] <- 0

        s2[l] <- sumaConPesos(estim.ini[l,],x,y,cw, w,uu)
    }

    s0  <- order(s2)
    beta0ini <- estim.ini[s0[1],]

    list(beta0ini,kk)
}
#####################################################################
betaExacto <- function(sbmn,x,y)
{
    ## to each subsample assign the maximum likelihood estimator and
    ## fixing the case mle has NA components
    nmue <- nrow(x)
    p <- ncol(x)
    betaexacto <- glm(y[sbmn]~x[sbmn,]-1,family = poisson())$coefficients
    sinNas <- na.exclude(betaexacto)
    long <- length(sinNas)
    lugaresNas <- na.action(sinNas)[1:(p-long)]
    betaexactoSinNas <- betaexacto
    betaexactoSinNas[lugaresNas] <- 0
    betaexactoSinNas
}



###################################################################################
glmRobMT <- function(x,y,cw = 2.1,iw = 0,nsubm = 500)
{
    ## MAINFUNCTION  Computes the MT or WMT estimator  for Poisson regression  with intercept starting from the estimator computed in the function
    ##beta0IniC.
    ##INPUT
    ##x  design   matrix with nrows and p columns.
    ##y respone  vector of lengtg n
    ##cw tuning constant. Default value 2.1
    ##iweigths indicator for weights penalizing high leverage points, iweights=1 indicates to use weights iweights=0
    ##indicate notto use way. Default value is iw=0, Our simulation study suggests not to use weights.
    ##nsubm Number of subsamples. Default calue nsubm=500
    ##OUTPUT
    ##$initial is the inital estimate (first component is the intercept)
    ##$final is the final estimate (first component is the intercept)
    ##$nsamples is the number of well  conditioned  subsamples
    ##REQUIRED PACKAGES: tools, rrcov
    uu <- spcreation(cw)
    n <- nrow(x)
    x1 <- cbind(rep(1,n),x)
    w <- weights(x1,iw)
    out <- beta0IniCP(x1,y,cw,w,uu,nsubm)
    beta0  <- out[[1]]
    nsubm  <- out[[2]]

    estim2 <- optim(beta0,sumaConPesos,method = "L-BFGS-B",x = x1,y = y,cw = cw, w = w,uu = uu)

    mtest <- estim2$par

    out  <- list("initial" = beta0, "final" = mtest,"nsubsamples" = nsubm)
    out
}



