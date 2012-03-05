lmrob.lar <- function(x, y, tol=1e-6)
{
  ## LAR : Least Absolute Residuals -- i.e. L_1  M-estimate
  ## this function is identical to lmRob.lar of the robust package

  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  stopifnot(p > 0, n >= p, length(y) == n, length(tol) == 1, is.numeric(tol))
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  bet0 <- 0.773372647623  ## bet0 = pnorm(0.75)
  tmpn <- double(n)
  tmpp <- double(p)

  z1 <- .Fortran(rllarsbi, ##-> ../src/rllarsbi.f
                 x,
                 y,
                 as.integer(n),
                 as.integer(p),
                 as.integer(n),
                 as.integer(n),
                 as.double(tol),
                 NIT=integer(1),
                 K=integer(1),
                 KODE=integer(1),
                 SIGMA=double(1),
                 THETA=tmpn,
                 RS=tmpn,
                 SC1=tmpn,
                 SC2=tmpp,
                 SC3=tmpp,
                 SC4=tmpp,
                 BET0=as.double(bet0),
                 PACKAGE = "robustbase")[c("THETA","SIGMA","RS","NIT")]
  names(z1) <- c("coef", "scale","resid","iter")
  ##           c("THETA","SIGMA", "RS",  "NIT")
  length(z1$coef) <- p
  z1
  ##list(coef=z1$THETA[1:p], scale=z1$SIGMA, resid=z1$RS)
}

lmrob.split <- function(x, y, control, mf) {
    mt <- attr(mf, "terms")
    p <- ncol(x)
    
    ## --- split categorical and interactions of categorical vars.
    ##     from continuous variables
    factors <- attr(mt, "factors")
    factor.idx <- attr(mt, "dataClasses") == "factor"
    if (!any(factor.idx)) ## There are no factors
        return(list(x1.idx=rep(FALSE, p), x2=x))
    ## --- include interactions cat * cont in x1:
    ## factor.asgn <- which(factor.idx %*% factors > 0)
    ## --- include also continuous variables that interact with factors in x1:
    ##     make sure to include interactions of continuous variables
    ##     interacting with categorical variables, too
    factor.asgn <- numeric(0)
    factors.cat <- factors[, factor.idx %*% factors > 0,drop=FALSE]
    factors.cat[factors.cat > 0] <- 1 ## fix triple+ interactions
    for (i in 1:ncol(factors)) {
        comp <- factors[,i] == 1
        ## if any of the components is a factor: include in x1 and continue
        if (any(factor.idx[comp])) {
            factor.asgn <- c(factor.asgn, i)
        } else {
            ## if there is an interaction of this term with a categorical var.
            ## include in x1 and continue
            if (any(colSums(factors.cat[comp,,drop=FALSE]) >= sum(comp)))
                factor.asgn <- c(factor.asgn, i)
        }
    }
    ## --- do not include interactions cat * cont in x1:
    ## factor.asgn <- which(factor.idx %*% factors & !(!factor.idx) %*% factors)
    x1.idx <- attr(x, "assign") %in% c(0, factor.asgn) ## also include intercept
    names(x1.idx) <- colnames(x)
    
    ## x1: factors and interactions of / with factors
    ## x2: continuous variables (incl. intercept)
    x1 <- x[, x1.idx, drop=FALSE]
    x2 <- x[, !x1.idx, drop=FALSE]
    
    p1 <- ncol(x1)
    p2 <- ncol(x2)

    if (p2 == 0)
        return(list(x1=x1, x1.idx=x1.idx))
    
    ## --- initial estimates: lar / L1
    tmp <- lmrob.lar(x1, y, control$rel.tol)
    y.tilde <- tmp$resid
    t1 <- tmp$coef
    x2.tilde <- x2
    T2 <- matrix(0, nrow=p1, ncol=p2)
    
    for(i in 1:p2) {
        tmp <- lmrob.lar(x1, x2[,i], control$rel.tol)
        x2.tilde[,i] <- tmp$resid
        T2[,i] <- tmp$coef
    }
    list(x1=x1, x1.idx=x1.idx, x2=x2, x2.tilde=x2.tilde,
         y.tilde=y.tilde, t1=t1, T2=T2)
}
