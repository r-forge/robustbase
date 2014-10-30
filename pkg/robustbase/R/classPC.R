##' @title Simple Matrix Rank ====> ../man/rankMM.Rd
rankMM <- function(A, tol = NULL, sv = svd(A,0,0)$d) {
    d <- dim(A)
    stopifnot(length(d)==2, length(sv)==min(d), diff(sv) <= 0) # must be sorted decreasingly
    if(is.null(tol))
	tol <- max(d) * .Machine$double.eps * abs(sv[1])
    else
	stopifnot(is.numeric(tol), tol >= 0)
    sum(sv >= tol)
}

##' Flip the signs of the loadings
##'  - comment from Stephan Milborrow
.signflip <- function(loadings) {
    apply(loadings, 2L,
	  function(x) if(x[which.max(abs(x))] < 0) -x else x)
}

##' @title Classical Principal Components ... ==> ../man/classPC.Rd
classPC <- function(x, scale=FALSE, center=TRUE,
		    signflip=TRUE, via.svd = n > p, scores=FALSE)
{
    if(!is.numeric(x) || !is.matrix(x))
	stop("'x' must be a numeric matrix")
    else if((n <- nrow(x)) <= 1)
	stop("The sample size must be greater than 1 for svd")
    p <- ncol(x)
    x <- scale(x, center=center, scale=scale)
    ##	 -----
    if(isTRUE(scale))
	scale <- attr(x, "scaled:scale")
    if(isTRUE(center))
	center <- attr(x, "scaled:center")

    if(via.svd) {
	svd <- svd(x/sqrt(n-1), nu=0)
	rank <- rankMM(x, sv=svd$d)
	loadings <- svd$v[,1:rank]
	eigenvalues <- (svd$d[1:rank])^2 ## FIXME: here .^2; later sqrt(.)
    } else { ## n <= p; was "kernelEVD"
	e <- eigen(tcrossprod(x)/(n-1), symmetric=TRUE)
	tolerance <- n * max(e$values) * .Machine$double.eps
	rank <- sum(e$values > tolerance)
	ii <- seq_len(rank)
	eigenvalues <- e$values[ii]
	## MM{FIXME (efficiency)}:
	loadings <- t((x/sqrt(n-1))) %*% e$vectors[,1:rank] %*% diag(1/sqrt(eigenvalues))
    }

    ## VT::15.06.2010 - signflip: flip the sign of the loadings
    if(signflip)
	loadings <- .signflip(loadings)

    list(rank=rank, eigenvalues=eigenvalues, loadings=loadings,
	 scores = if(scores) x %*% loadings,
	 center=center, scale=scale)
}
