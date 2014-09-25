rlas <- function(y, b=0.2, mfn = function(n) 0.1*n^(-0.25),
                 nstart=30, scon=NULL)
{
    ## Initialize:
    y0    <- y[1:nstart]
    alpha <- median(y0)
    s     <- if(is.null(scon)) mean(abs(y0-alpha)) else scon
    mu    <- mfn(nstart)/s
    eps   <- s/nstart^b
    kount <- sum(abs(alpha-y0) < eps)
    g     <- kount/(eps*nstart)
    ny    <- length(y)
    n     <- nstart+1
    locn  <- numeric(ny)
    locn[1:(nstart-1)] <- NA
    locn[nstart] <- alpha
    scale  <- numeric(ny)
    scale[1:(nstart-1)] <- NA
    scale[nstart] <- s

    ## Calculate recursively, i.e., update (s, g, alpha)
    ## and store all (alpha, scale) in (locn,scale)
    while(n <= ny) {
        d.n <- y[n] - alpha
        s     <- if(is.null(scon)) ((n-1)*s + abs(d.n))/n else scon
        mu    <- mfn(n)/s
        eye   <- if(abs(d.n) < s/n^b) 1 else 0
        g     <- (1-1/n)*g + n^(b-1)*eye/s
        a     <- max(mu, g)
        alpha <- alpha + sign(d.n)/(a*n)
        locn[n]  <- alpha
        scale[n] <- s
        n <- n+1
    }
    list(locn=locn, scale=scale)
}
