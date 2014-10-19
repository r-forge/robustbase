##  .detmcd()
##
##  .detmcd computes the MCD estimator of a multivariate data set
##  in a deterministic way.
##  This estimator is given by the subset of h observations with smallest
##  covariance determinant.  The MCD location estimate is then
##  the mean of those h points, and the MCD scatter estimate is
##  their covariance matrix.  The default value of h is roughly
##  0.75n (where n is the total number of observations), but the
##  user may choose each value between n/2 and n. Based on the
##  raw estimates, weights are assigned to the observations such
##  that outliers get zero weight. The reweighted MCD estimator
##  is then given by the mean and covariance matrix of the cases
##  with non-zero weight.
##
##  To compute the MCD estimator, six initial robust h-subsets are
##  constructed based on robust transformations of variables or robust and
##  fast-to-compute estimators of multivariate location and shape. Then
##  C-steps are applied on these h-subsets until convergence. Note that the
##  resulting algorithm is not fully affine equivariant, but it is often
##  faster than the FAST-MCD algorithm which is affine equivariant
##  (see covMcd()).
##  Note that this function can not handle exact fit situations: if the
##  raw covariance matrix is singular, the program is stopped. In that
##  case, it is recommended to apply the covMcd() function.
##
## Reference:
##  Hubert, M., Rousseeuw, P.J. and Verdonck, T. (2012),
##  "A deterministic algorithm for robust location and scatter", Journal of
##  Computational and Graphical Statistics, in press.
##
##  The MCD method is intended for continuous variables, and assumes that
##  the number of observations n is at least 5 times the number of variables p.
##  If p is too large relative to n, it would be better to first reduce
##  p by variable selection or robust principal components (see the functions
##  robust principal components in package 'rrcov').
##
## Input argument:
##  x    - a numerical matrix. The columns represent variables, and rows represent observations.
##  h    - The quantile of observations whose covariance determinant will
##          be minimized.  Any value between n/2 and n may be specified.
##          The default value is 0.5*n.
## scales-  function to estimate the scale. Default is Qn, the user
##          could choose mad or scaleTau2, or anything else
##
## hsets.init: If one gives here already a matrix with for each row an
##             ordering of the observations (first the one with smallest statistical
##             distance), then the initial shape estimates are not calculated.
##             Default value = NULL.

.detmcd <- function(x,
                    h,
                    hsets.init=NULL,
                    scales,
                    trace=as.integer(trace))
{
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    kmini <- 5      # number of sub-data sets(if we use them some day)
                    # for now we use it as number of rows in the returned
                    # matrix 'coeff' for exact fit (also not used currently).

    if(missing(scales) | is.null(scales))
    {
        scales <- if(n <= 5000) Qn else scaleTau2
    }

    cutoff = qchisq(0.975, p)
    chimed = qchisq(0.5,   p)

    ## Center and scale the data
    z <- doScale(x, center=median, scale=scales)
    z.center <- z$center
    z.scale <- z$scale
    z <- z$x

    if(is.null(hsets.init))
        hsets.init <- getInitSubsets(z, scaled=TRUE, h=h, scales=scales)

    ## Construction of h-subset: take the first h observations
    hsets <- hsets.init[, 1:h, drop=FALSE]       ## select the first h ranked observations
    nsets <- NROW(hsets)

    ## Some initializations.
    raw.wt <- rep(NaN, n)
    raw.rd <- rep(NaN, n)
    rew.rd <-rep(NaN, n)
    rew.mahalanobis <- rep(NaN, n)
    rew.flag <- rep(NaN, n)

    rew.Hsubsets.Hopt <- NULL
    rew.Hsubsets.i <- NULL
    rew.Hsubsets.csteps <- rep(NA, nsets)

    csteps  <- 200                          # 100 csteps ?
    prevdet <- 0                            # ?
    bestobj <- Inf                          # inf in R ?
    teller  <- rep(0, n+1)

    for(i in 1:nsets)
    {
        for(j in 1:csteps)
        {
            if(j == 1)  {
                obs_in_set <- hsets[i,]     # start with the i-th initial set
            } else
            {
                score <- (z - repmat(svd$center, n, 1)) %*% svd$loadings
                mah <- mahalanobis(score, matrix(0, p, 1), diag(svd$eigenvalues))
                ord <- order(mah);
                obs_in_set <- ord[1:h]
            }
            svd <- classSVD(z[obs_in_set,])      # [P,T,L,r,centerX,meanvct] = classSVD(data(obs_in_set,:));
            obj <- prod(svd$eigenvalues)

            if(svd$rank < p)
                stop('More than h of the observations lie on a hyperplane.')
            if(j >= 2 && obj == prevdet)
                break
            prevdet <- obj
        }

        if(obj < bestobj)
        {
            ## bestset           : the best subset for the whole data.
            ## bestobj           : objective value for this set.
            ## initmean, initcov : resp. the mean and covariance matrix of this set

            bestset <- obs_in_set
            bestobj <- obj
            initmean <- svd$center
            initcov <- svd$loadings %*% diag(svd$eigenvalues) %*% t(svd$loadings)
            raw.initcov <- initcov
            rew.Hsubsets.Hopt <- bestset
            rew.Hsubsets.i <- i         # to determine which subset gives best results.
        }

        rew.Hsubsets.csteps[i] <- j     # how many csteps necessary to converge.
    }

    svd <- classSVD(z[bestset, ])    # [P,T,L,r,centerX,meanvct] = classSVD(data(bestset,:));
    mah <- mahalanobis((z - repmat(svd$center, n, 1)) %*% svd$loadings, matrix(0, p, 1), diag(svd$eigenvalues))
    sortmah <- sort(mah)

##    factor <- sortmah[h]/qchisq(h/n, p)
##    raw.cov <- factor*initcov

    raw.cov <- initcov
    ## We express the results in the original units.
    raw.cov <- raw.cov * repmat(z.scale, p, 1) * t(repmat(z.scale, p, 1))
    raw.center <- initmean * z.scale + z.center
    raw.objective <- bestobj * prod(z.scale)^2
    raw.mah <- mahalanobis(x, raw.center, raw.cov, tol=1E-14)
    medi2 <- median(raw.mah)

    ans <- list(initcovariance=raw.cov,
                initmean=raw.center,
                best=sort(bestset),
                mcdestimate=raw.objective,  # determinant (goes to crit)
                weights=NULL,               # goes to raw.weights
                exactfit=0,
                coeff=matrix(rep(0, kmini*p), nrow=kmini),
                kount=0)
}

doScale <- function (x, center=mean, scale=sd)
{
    stopifnot(class(x) == "matrix");
	n <- nrow(x)
	p <- ncol(x)

	m <- array(0, p)
	if(is.character(center))
		center = eval(parse(text=center))

    if(is.function(center))
	{
		m <- apply(x, 2, center)
		x <- x - matrix(1, nrow = n) %*% m
	}else if(length(center) == p & is.numeric(center))
	{
		m <- center
		x <- x - matrix(1, nrow = n) %*% m
	}else if(!is.null(center))
		warning ("Unknown center specification! Centering will be omitted.\n")

	s <- array(1, p)
	if(is.character(scale))
		scale <- eval(parse(text=scale))

	if(is.function(scale))
	{
		s <- apply(x, 2, scale)
		x <- x / matrix(1, nrow = n) %*% s
	}else if(length(scale) == p & is.numeric(scale))
	{
		s <- scale
		x <- x / matrix(1, nrow = n) %*% s
	}else if(!is.null(scale))
		warning ("Unknown scale specification! Scaling will be omitted.\n")

	return (list(x=x, center=m, scale=s))
}

classSVD <- function(x, scale=FALSE, signflip=TRUE){

    ## Flip the signs of the loadings
    ##  - comment from Stephan Milborrow
    ##
    .signflip <- function(loadings)
    {
        if(!is.matrix(loadings))
            loadings <- as.matrix(loadings)
        apply(loadings, 2, function(x) if(x[which.max(abs(x))] < 0) -x else x)
    }

    if(!is.numeric(x) || !is.matrix(x))
        stop("'x' must be a numeric matrix")
    else if(nrow(x) <= 1)
        stop("The sample size must be greater than 1 for svd")

    n <- nrow(x)
    p <- ncol(x)

    center <- apply(x, 2, mean)
    x <- scale(x, center=TRUE, scale=scale)
    if(scale)
        scale <- attr(x, "scaled:scale")

    svd <- svd(x/sqrt(n-1))

    rank <- rankMM(x, sv=svd$d)
    eigenvalues <- (svd$d[1:rank])^2
    loadings <- svd$v[,1:rank]

    ## VT::15.06.2010 - signflip: flip the sign of the loadings
    if(!is.matrix(loadings))
        loadings <- data.matrix(loadings)
    if(signflip)
        loadings <- .signflip(loadings)
    scores <- x %*% loadings

    list(loadings=loadings,
         scores=scores,
         eigenvalues=eigenvalues,
         rank=rank,
         x=x,
         center=center,
         scale=scale)
}

##  <Matlab>
##      a=[1 2 ; 3 4];
##      repmat(a,2,3)
##
##  <R>
##      a <- matrix(1:4,2,byrow=T)
##      repmat(a,2,3)

repmat <- function(A, n, p) {

    if(is.vector(A))    # we need a column matrix, not a vector, speaking in R terms
        A <- t(A)
    kronecker(matrix(1,n,p), A)
}

## Purpose: rank of a matrix ``as Matlab''
## ----------------------------------------------------------------------
## Arguments:  A: a numerical matrix, maybe non-square
##           tol: numerical tolerance (compared to singular values)
##            sv: vector of non-increasing singular values of A
##                (pass as argument if already known)
## ----------------------------------------------------------------------
## Author: Martin Maechler, Date:  7 Apr 2007, 16:16
rankMM <- function(A, tol = NULL, sv = svd(A,0,0)$d) {
    d <- dim(A)
    stopifnot(length(d)==2, length(sv)==min(d), diff(sv) <= 0)   # must be sorted decreasingly
    if(is.null(tol))
        tol <- max(d) * .Machine$double.eps * abs(sv[1])
    else
        stopifnot(is.numeric(tol), tol >= 0)
    sum(sv >= tol)
}

##  Compute six initial (robust) estimators of location scale,
##  for each of them compute the distances and take the h0 = n/2
##  observations with smallest dik. Then compute the statistical
##  distances based on these h0 observations. Return the indexes
##  of the observations sorted in accending order.
##
getInitSubsets <- function(x, scaled=TRUE, h, scales)
{
    ##
    ## As the considered initial estimators Sk may have very
    ##  inaccurate eigenvalues, we try to 'improve' them by applying
    ##  transformation similar to that used in the OGK algorithm.
    ##  After that compute the corresponding distances, ordther them and
    ##  return the indexies
    ##
    initset <- function(data, scales, P)
    {
        n <- NROW(data)
        p <- NCOL(data)
        lambda <- doScale(data %*% P, center=median, scale=scales)$scale
        sqrtcov <- P %*% diag(lambda) %*% t(P)
        sqrtinvcov <- P %*% diag(1/lambda) %*% t(P)
        estloc <- apply(data %*% sqrtinvcov,2, median) %*% sqrtcov
        centeredx <- (data - repmat(estloc, n, 1)) %*% P
        order(mahalanobis(centeredx, matrix(0, p, 1), diag(lambda)^2))
    }

    ##
    ##  Compute the raw OGK estimator. For m(.) and s(.) (robust
    ##  univariate estimators of location and scale) use the median
    ##  and Qn for reasons of simplicity (no choice of tuning parameters)
    ##  and to be consistent with the other components of DetMCD.
    ##
    ogkscatter <- function(Y, scales)
    {
        if(!is.matrix(Y))
            Y <- as.matrix(Y)

        n <- nrow(Y)
        p <- ncol(Y)
        U <- diag(p)

        for(i in 1:p)
        {
            if( i > 1)
            {
                sYi <- Y[,i]
                for(j in 1:(i-1))
                {
                    sYj <- Y[,j]
                    ssy <- scales(sYi + sYj)
                    sdy <- scales(sYi - sYj)
                    U[i,j] <- (ssy^2 - sdy^2) / 4
                }
            }
        }

        U <- lower.tri(U) * U + t(U)    #    U <- tril(U, -1) + t(U)
        ee <- eigen(U)
        P <- ee$vectors

        Z <- Y %*% t(P)
        sigz <- apply(Z, 2, scales)
        lambda <- diag(sigz^2)

        list(P=P, lambda=lambda)
    }

    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]

    if(missing(scales) | is.null(scales))
    {
        scales <- if(n <= 5000) Qn else scaleTau2
    }

    ## If the data was not scaled (scaled=FALSE), center and scale using
    ##  the median and the provided function 'scales'. If scales is missing
    ##  or is NULL, use Qn, for smaller data sets (n <= 5000) and
    ##  tau-scale of Yohai and Zamar (1988) otherwise (implemented in
    ##  the function scaleTau2() in 'robustbase').
    z <- if(!scaled)
        {
            ## Center and scale the data
            z <- doScale(x, center=median, scale=scales)
            z.center <- z$center
            z.scale <- z$scale
            z <- z$x
        } else
            x

    hsets.init=NULL

    ## Determine 6 initial estimates (ordering of obs)
    ## 1. Hyperbolic tangent of standardized data
    y1 <- tanh(z)
    R1 <- cor(y1)
    P <- eigen(R1)$vectors
    ind <- initset(z, scales=scales, P=P)
    hsets.init <- rbind(hsets.init, ind)

    ## 2. Spearmann correlation matrix
    R2 <- cor(z, method="spearman")
    P <- eigen(R2)$vectors
    ind <- initset(z, scales=scales, P=P)
    hsets.init <- rbind(hsets.init, ind)

    ## 3. Tukey normal scores
    y2 <- z
    for(j in 1:p)
        y2[,j] <- rank(z[, j])
    y3 <- qnorm((y2-1/3)/(n+1/3))
    R3 <- cor(y3, use = "complete.obs")
    P <- eigen(R3)$vectors
    ind <- initset(z, scales=scales, P=P)
    hsets.init <- rbind(hsets.init, ind)

    ## 4. Spatial sign covariance matrix
    znorm <- as.matrix(sqrt(apply(z^2, 1, sum)))
    ii <- znorm > .Machine$double.eps
    zznorm <- z
    zznorm[ii,] <- z[ii, ] / repmat(znorm[ii,, drop=FALSE], 1, p)
    SCM <- (t(zznorm) %*% zznorm) / (n-1)
    P <- eigen(SCM)$vectors
    ind <- initset(z, scales=scales, P=P)
    hsets.init <- rbind(hsets.init, ind)

    ## 5. BACON
    ind5 <- order(znorm)
    half <- ceiling(n/2);
    Hinit <- ind5[1:half]
    covx <- cov(z[Hinit,])
    P <- eigen(covx)$vectors
    ind <- initset(z, scales=scales, P=P)
    hsets.init <- rbind(hsets.init, ind)

    ## 6. Raw OGK estimate for scatter
    P <- ogkscatter(z, scales)$P
    ind <- initset(z, scales=scales, P=P)
    hsets.init <- rbind(hsets.init, ind)

    hsets.init <- as.matrix(hsets.init)

    hsets <- hsets.init[, 1:h, drop=FALSE]       ## select the first h ranked observations
    nsets <- NROW(hsets)

    for(k in 1:nsets)
    {
        xk <- z[hsets[k,] ,]
        svd <- classSVD(xk)         # [P,T,L,r,centerX,meanvct] = classSVD(xk)
        if(svd$rank < p)
            stop('DetMCD.m: More than half of the observations lie on a hyperplane.')

        score <- (z - repmat(svd$center, n, 1)) %*% svd$loadings
        ord <- order(mahalanobis(score, matrix(0, p, 1), diag(svd$eigenvalues)))
        hsets.init[k,] <- ord
    }

    hsets.init
}
