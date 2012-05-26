### test subsample
### LU decomposition and singular subsamples handling
require(robustbase)
source(system.file("xtraR/subsample-fns.R", package = "robustbase", mustWork=TRUE))

gotMatrix <- require(Matrix) ## to test lu decomposition

## R version of LU decomposition in subsample() in lmrob.c
## Modified from Golub G. H., Van Loan, C. F., Matrix Computations, 3rd edition
LU.gaxpy <- function(A, pivot=TRUE) {
    A <- as.matrix(A)
    ## precondition:
    cf0 <- max(abs(A))
    A <- A / cf0
    stopifnot(ncol(A) >= nrow(A))
    m <- nrow(A)
    n <- ncol(A)
    v <- double(m) ## work matrix
    ## these matrices will contain the results
    L <- diag(m)
    U <- matrix(0, m, m)
    p <- integer(m-1) ## pivots
    idc <- 1L:n ## which columns of A are used
    idr <- 1L:m ## how to rows of A are permuted
    lenidc <- length(idc)
    for(j in 1L:m) {
        sing <- TRUE
        while(sing) {
            if (j == 1L) {
                v[j:m] <- A[idr[j:m], idc[j]]
            } else {
                rows <- 1L:(j-1L)
                z <- forwardsolve(L[rows, rows, drop=FALSE], A[idr[rows], idc[j]])
                ## cat("Step", j, "z=", sapply(z, function(x) sprintf("%.15f", x)), "\n");
                U[rows, j] <- z
                v[j:m] <- A[idr[j:m], idc[j]] - L[j:m, rows, drop=FALSE] %*% z
                ## cat("v=", v, "\n");
            }
            if (j < m) {
                mu <- j
                mu <- if (pivot) which.max(abs(v[j:m])) + j - 1L else j
                ## debug possumDiv example
                ## cat(sprintf("R-Step: %i: ", j), round(abs(v[j:m]), 6), "\n"); 
                ## cat(mu, v[mu], "\n")
                if (abs(v[mu]) >= 1e-7) { ## singular: can stop here already
                    p[j] <- mu
                    if (pivot) {
                        tmp <- v[j]; v[j] <- v[mu]; v[mu] <- tmp
                        tmp <- idr[j]; idr[j] <- idr[mu]; idr[mu] <- tmp
                    }
                    L[(j+1L):m, j] <- v[(j+1L):m]/v[j]
                    if (pivot & j > 1) {
                        tmp <- L[j, rows]; L[j, rows] <- L[mu, rows]; L[mu, rows] <- tmp
                    }
                }
            }
            U[j, j] <- v[j]
            if (abs(v[j]) < 1e-7) {
                ##warning("singularity detected in step ", j, " candidate: ", idc[j])
                idc[j] <- idc[lenidc]
                lenidc <- lenidc - 1
                if (lenidc < j)
                    break
            } else sing <- FALSE
        }
    }
    list(L = L, U = U * cf0, p = p, idc = idc[1L:m], singular = sing)
}

Rsubsample <- function(x, y, mts=0) {
    x <- as.matrix(x)
    n <- length(y)
    p <- ncol(x)

    .C(robustbase:::R_subsample,
       x=as.double(x),
       y=as.double(y),
       n=n,
       m=p,
       beta=double(p),
       ind_space=integer(n),
       idc=integer(n),
       idr=integer(n),
       lu=double(p^2),
       v=double(p),
       pivot=integer(p-1),
       status=integer(1),
       sample=FALSE,
       mts=as.integer(mts),
       ss=as.integer(mts == 0))
}

subsample <- function(x, y=rnorm(nrow(x))) {
    x <- as.matrix(x)
    n <- length(y)
    p <- ncol(x)

    z <- Rsubsample(x, y)
    ## convert idc, idr and p to 1-based indexing
    idr <- z$idr + 1
    idc <- z$idc[1:p] + 1
    pivot <- z$p + 1
    ## get L and U
    L <- U <- LU <- matrix(z$lu, p, p)
    L[upper.tri(L, diag=TRUE)] <- 0
    diag(L) <- 1
    U[lower.tri(U, diag=FALSE)] <- 0

    cmp <- LU.gaxpy(t(x))
    if (!isTRUE(all.equal(cmp$p, pivot))) {
        print(cmp$p)
        print(pivot)
        print(which(cmp$p != pivot))
    }
    stopifnot(all.equal(cmp$L, L),
              all.equal(cmp$U, U),
              all.equal(cmp$p, pivot),
              all.equal(cmp$idc, idc))

    xsub <- x[idc, ]
    ysub <- y[idc]

    ## compare with Matrix result
    if (gotMatrix & !cmp$singular) {
        cf0 <- max(abs(x))
        tmp <- lu(t(xsub / cf0))
        idx <- upper.tri(xsub, diag=TRUE)
        tmp@x[idx] <- tmp@x[idx] * cf0
        stopifnot(all.equal(tmp@x, z$lu))
    }

    ## test solved parameter
    if (!cmp$singular) {
        stopifnot(all.equal(z$beta, unname(solve(xsub, ysub))))
    }

    invisible(z)
}

A <- matrix(c(0.001, 1, 1, 2), 2)
set.seed(11)
str(sa <- subsample(A))

A <- matrix(c(3, 2, 6, 17, 4, 18, 10, -2, 12), 3)
subsample(A)

## test some random matrix
set.seed(1002)
A <- matrix(rnorm(100), 10)
subsample(A)

## test singular matrix handling
A <- matrix(c(1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1), 4, byrow=TRUE)
subsample(A)


## test subsample with mts > 0
data <- data.frame(y = rnorm(9), expand.grid(A = letters[1:3], B = letters[1:3]))
x <- model.matrix(y ~ ., data)
y <- data$y
## this should produce a warning and return status == 2
z <- Rsubsample(x, y, mts=2)
stopifnot(z$status == 2)

## test real data example
data(possumDiv)## 151 * 9; the last two variables are factors
with(possumDiv, table(eucalyptus, aspect))

mf <- model.frame(Diversity ~ .^2, possumDiv)
X <- model.matrix(mf, possumDiv)
y <- model.response(mf)
stopifnot(qr(X)$rank == ncol(X))

## test pre-conditioning
## this used to fail: different pivots in step 37
str(s1 <- subsample(X, y))
s2 <- subsample(X / max(abs(X)), y / max(abs(X)))
s3 <- subsample(X * 2^-50, y * 2^-50)
## all components *BUT*  x, y, lu :
nm <- names(s1); nm <- nm[is.na(match(nm, c("x","y","lu")))]
stopifnot(all.equal(s1[nm], s2[nm], tol=1e-10),
	  all.equal(s1[nm], s3[nm], tol=1e-10))

## test subsampling
testSubSampling <- function(X, y) {
    lX <- X[sample(nrow(X)), ]
    ## C version
    zc <- Rsubsample(lX, y)
    ## R version
    zR <- LU.gaxpy(t(lX))
    if (as.logical(zc$status)) {
        ## singularity in C detected
        if (!zR$singular)
            stop("singularity in C but not in R")
    } else {
        ## no singularity detected
        if (zR$singular)
            stop("singularity in R but not in C")
    }
    zR$singular
}

set.seed(10)
nsing <- sum(replicate(200, testSubSampling(X, y)))
stopifnot(nsing == 0)

## test example with many categorical predictors
set.seed(10)
r1 <- lmrob(Diversity ~ .^2 , data = possumDiv, cov="none")
## lmrob.S used to fail for this seed:
set.seed(108)
lmrob(Diversity ~ .^2 , data = possumDiv, cov="none") #, trace=4)

## investigate problematic subsample:
idc <- 1 + c(140, 60, 12, 13, 89, 90, 118, 80, 17, 134, 59, 94, 36,
         43, 46, 93, 107, 62, 57, 116, 11, 45, 35, 38, 120, 34, 29,
         33, 147, 105, 115, 92, 61, 91, 104, 141, 138, 129, 130, 84,
         119, 132, 6, 135, 112, 16, 67, 41, 102, 76, 111, 82, 148, 24,
         131, 10, 96, 0, 87, 21, 127, 56, 124)

rc <- lm(Diversity ~ .^2 , data = possumDiv, subset = idc)

X <- model.matrix(rc)
y <- possumDiv$Diversity[idc]
subsample(X, y)

lu <- LU.gaxpy(t(X))
stopifnot(lu$sing)
zc <- Rsubsample(X, y)
stopifnot(zc$status > 0)
## column 52 is linearly dependent and should have been discarded
## qr(t(X))$pivot

image(as(round(zc$lu - (lu$L + lu$U - diag(nrow(lu$U))), 10), "Matrix"))
image(as(sign(zc$lu) - sign(lu$L + lu$U - diag(nrow(lu$U))), "Matrix"))
