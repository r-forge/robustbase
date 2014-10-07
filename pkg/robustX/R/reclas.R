reclas <- function(y, b=0.2, mfn = function(n) 0.1*n^(-0.25),
                   nstart=30, m0 = median(y0),
                   scon=NULL, updateScale = is.null(scon))
{
    ## Initialize:
    y0 <- y[1:nstart]
    force(updateScale)
    alpha <- m0
    s <- if(is.numeric(scon) && length(scon) == 1 && scon > 0) scon
         else if(is.function(scon)) scon(y0, m0)
         else if(is.null(scon)) mean(abs(y0-alpha))
         else stop("Invalid 'scon': must be NULL, a function or a positive number")
    ny    <- length(y)
    locn  <- numeric(ny); locn [1:(nstart-1)] <- NA; locn [nstart] <- alpha
    if(updateScale) {
	scale <- numeric(ny); scale[1:(nstart-1)] <- NA; scale[nstart] <- s
    } else { # really not storing full length constant scale[] vector
	scale <- s
    }
    ## mu    <- mfn(nstart)/s
    eps   <- s/nstart^b
    kount <- sum(abs(alpha-y0) < eps)
    g     <- kount/(eps*nstart)
    n     <- nstart+1
    ## Calculate recursively, i.e., update (s, g, alpha)
    ## and store all (alpha, scale) in (locn,scale)
    while(n <= ny) {
        d.n <- y[n] - alpha
        if(updateScale) s <- ((n-1)*s + abs(d.n))/n
        s.n <- s / n^b
        mu    <- mfn(n)/s
        g     <- ((n-1)*g + if(abs(d.n) < s.n) 1/s.n else 0)/n
        a     <- max(mu, g)
        alpha <- alpha + sign(d.n)/(a*n)
        locn[n] <- alpha
        if(updateScale) scale[n] <- s
        n <- n+1
    }
    structure(list(locn=locn, scale=scale, updateScale=updateScale,
                   nstart=nstart, call = match.call()), # <- so methods can "display" the call
              class = "reclas")
}

plot.reclas <-
    function(x, M = tail(x$locn, 1), ylim = NULL,
             s.y = tail(x$scale, 1), se = TRUE, col.se = adjustcolor("skyblue4", 3/4),
             ylab = "locn", main = deparse(x$call)[1], ...)
{
   ## Default: use approximately *same* y-scale for the plots
   if(is.null(ylim)) ylim <- M + c(-.5,.5) * s.y
   L <- x$locn
   plot(L, type="l", ylim=ylim, ylab=ylab, main=main, ...)
   if(se) {
     n <- seq_along(L); sd <- x$scale / sqrt(n)
     lines(n, L + 2*sd, col=col.se)
     lines(n, L - 2*sd, col=col.se)
   }
   abline(h=M, col=adjustcolor("red",0.6), lty=2)
}
