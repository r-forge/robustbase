### TODOs / Questions (Martin  <-->  Eduardo):
### -------------------------------------------

## (2) argument names:
## 5 matches for "@param .*tuning" in buffer  nlregrob.R
##      29:##' @param tuning.chi.scale \
##      30:##' @param tuning.chi.M     /  nlrob.MM
##     208:##' @param tuning.chi.scale \
##     209:##' @param tuning.chi.tau   /  nlrob.tau
##     356:##' @param tuning.chi          nlrob.CM


##' Compute an MM-estimator for nonlinear (constrained) regression.
##'
##' Copyright 2013, Eduardo L. T. Conceicao
##' Available under the GPL (>= 2)
##' @title MM-estimator for Nonlinear Regression
##' @param formula
##' @param data
##' @param pnames
##' @param lower
##' @param upper
##' @param tolerance
##' @param f0
##' @param estim
##' @param psi
##' @param tuning.chi.scale
##' @param tuning.chi.M
##' @param optim.control
##' @param ... optional arguments for optimization, e.g., \code{trace = TRUE}, passed to \code{\link{JDEoptim}(.)}.
##' @return
##' @references
##' Yohai, V.J. (1987)
##' High breakdown-point and high efficiency robust estimates for regression.
##' \emph{The Annals of Statistics} \bold{15}, 642--656.
##' @examples
##' nlrob.MM( density ~ Asym/(1 + exp(( xmid -log(conc) )/scal )),
##'           data = DNase[DNase$Run == 1, ],
##'           pnames = c("Asym", "xmid", "scal"),
##'           lower = 0, upper = 3,
##'           optim.control = list(trace = 1), trace = TRUE )
##' @author Eduardo L. T. Conceicao
nlrob.MM <- function(formula, data, pnames, lower, upper,
                     tolerance = 1e-5, f0 = NULL,
                     init = c("S", "lts"),
                     psi = c("bisquare", "lqq", "optimal", "hampel"),
                     tuning.chi.scale = NULL, tuning.chi.M = NULL,
                     optim.control = list(), ...)
{
    init <- match.arg(init)
    psi <- match.arg(psi)
    if (is.null(tuning.chi.scale))
        tuning.chi.scale <- switch(psi, bisquare = 1.56,
                                   lqq = c(0.402, 0.268, 1.5),
                                   optimal = 0.405,
                                   hampel = c(1.5, 3.5, 8)*0.212)
    if (is.null(tuning.chi.M))
        tuning.chi.M <- switch(psi, bisquare = 4.68, lqq = c(1.473, 0.982, 1.5),
                               optimal = 1.060, hampel = c(1.5, 3.5, 8)*0.9014)

    c1 <- tuning.chi.scale
    c2 <- tuning.chi.M

    ## Preliminary psi-specific checks / computations:
    switch(psi,
           "bisquare" = {
               if(!require("PolynomF"))
		   stop(gettextf(
		"You must install package 'PolynomF' before you can use %s",
		"nlrob.MM(*, psi=bisquare)"), domain=NA)
           },
           "lqq" = { # lqqMax = rho(Inf), used in rho.inv() *and* ..
               c12 <- c1[1]+c1[2]
               lqqMax <- (c1[1]*c1[3] - 2*c12)/(1-c1[3]) + c12})

    rho1 <- function(t) Mchi(t, c1, psi)
    rho2 <- function(t) Mchi(t, c2, psi)

    rho.inv <- switch(psi, bisquare = function(y) {
        x <- polynom()
        p <- 3*x^2 - 3*x^4 + x^6
        t <- solve(p, y)
        t <- Re( t[Im(t) == 0] )
        unique( t[t >= 0] ) * c1
    }, lqq = function(y) {
        uniroot( function(x) rho1(x) - y, lower = 0, upper = lqqMax )$root
    }, optimal = function(y) {
        ## Salibian-Barrera, Matias, Willems, Gert, and Zamar, Ruben (2008).
        ## The fast-tau estimator for regression.
        ## Journal of Computational and Graphical Statistics 17, 659-682.
        sqrt(y/1.38) * c1 * 3
    }, hampel = function(y) {
        C <- MrhoInf(c1, psi)
        a <- c1[1]; b <- c1[2]; r <- c1[3]
        if (y <= a/C)
            sqrt(2*C*y)
        else if (y <= (2*b - a)/C)
            0.5*a + C/a*y
        else r + sqrt( r^2 - ( (r - b)*(2*C/a*y + (b - a)) - b*r ) )
    })

    M_scale <- function(sigma, u) sum( rho1(u/sigma) )/nobs - 0.5

    objective.initial <-
        switch(init,
               "lts" = function(par) { ## 'h' is defined "global"ly
                   names(par) <- pnames
                   fit <- eval( formula[[3L]], c(data, par) )
                   sum(sort.int( (y - fit)^2, partial = h )[1:h])
               },
               "S" = function(par) {
                   names(par) <- pnames
                   fit <- eval( formula[[3L]], c(data, par) )
                   res <- y - fit
                   ## Rousseeuw, Peter J., and Leroy, Annick M. (1987).
                   ## Robust Regression & Outlier Detection.
                   ## John Wiley & Sons, New York, p. 137.
                   med_abs_res <- median(abs(res))
                   sigma <- uniroot( M_scale,
                                    lower = constant[1L] * med_abs_res,
                                    upper = constant[2L] * med_abs_res,
                                    u = res )$root
                   sigma
               })

    objective.M <- function(par, sigma) {
        names(par) <- pnames
        fit <- eval( formula[[3L]], c(data, par) )
        sum(rho2( (y - fit)/sigma ))
    }

    formula <- as.formula(formula)
    dataName <- substitute(data)
    varNames <- all.vars(formula)
    obsNames <- rownames(data <- as.data.frame(data))
    data <- as.list(data)# to be used as such
    if (length(formula) == 2L) { ## as nls
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 0
    }

    if (!is.character(pnames) || any(is.na(match(pnames, varNames))))
        stop("'pnames' must be a character vector named with parameters in 'formula'")
    npar <- length(pnames)
    if (length(lower) == 1)
        lower <- rep(lower, npar)
    if (length(lower) != npar)
        stop(gettextf("lower must be either of length %d, or length 1", npar))
    if (length(upper) == 1)
        upper <- rep(upper, npar)
    if (length(upper) != npar)
        stop(gettextf("upper must be either of length %d, or length 1", npar))
    stopifnot(is.numeric(lower), is.numeric(upper), lower <= upper)
    y <- eval(formula[[2L]], data)
    nobs <- length(y)
    stopifnot(nobs >= npar)
    if (is.null(f0))
        f0 <- sum(scale(y, scale = FALSE)^2)
    stopifnot(is.numeric(f0), f0 > 0)
    constant <- c(
        switch(psi, bisquare = 1/c1,
               lqq = 1/lqqMax,
               optimal = 1/c1 * 1/3,
               hampel = 1/c1[3]),
        if(nobs %% 2) 2/rho.inv(2/(nobs+2)) else 1/rho.inv(1/(nobs+1)))
    switch(init, lts = h <- (nobs + npar + 1)%/%2)

    initial <- JDEoptim(lower, upper, objective.initial,
                        tol = tolerance, fnscale = f0, ...)
    names(initial$par) <- pnames
    res <- y - eval( formula[[3L]], c(data, initial$par) )

    med_abs_res <- median(abs(res))
    sigma <- uniroot( M_scale,
                      lower = constant[1L] * med_abs_res,
                      upper = constant[2L] * med_abs_res,
                      u = res )$root

    con <- list(fnscale = initial$value, parscale = initial$par)
    M <- optim( initial$par, objective.M, sigma = sigma,
                method = "L-BFGS-B", lower = lower, upper = upper,
                control = c(con, optim.control) )
    coef <- M$par
    names(coef) <- pnames
    status <-
	if (M$convergence == 0) "converged"
	else if (M$convergence == 1)
	    "maximum number of iterations reached without convergence"
	else M$message
    fit <- eval( formula[[3L]], c(data, coef) )
    names(fit) <- obsNames
    structure(list(call = match.call(), formula=formula, nobs=nobs,
                   coefficients = coef,
                   fitted.values = fit,
                   residuals = y - fit,
                   crit = M$value,
                   initial = initial,
                   Scale = sigma,
                   status = status, counts = M$counts, data = dataName),
              class = "nlrob")
} ## nlrob.MM


##' Compute an Tau-estimator for nonlinear (constrained) regression.
##'
##' Copyright 2013, Eduardo L. T. Conceicao
##' Available under the GPL (>= 2)
##' @title Tau-estimator for Nonlinear (Constrained) Regression
##' @param formula
##' @param data
##' @param pnames
##' @param lower
##' @param upper
##' @param tolerance
##' @param f0
##' @param psi
##' @param tuning.chi.scale
##' @param tuning.chi.tau
##' @param ... optional arguments for optimization, e.g., \code{trace = TRUE}, passed to \code{\link{JDEoptim}(.)}.
##' @return
##' @references
##' Yohai, V.J., and Zamar, R.H. (1988).
##' High breakdown-point estimates of regression by means of the minimization
##' of an efficient scale.
##' \emph{Journal of the American Statistical Association} \bold{83}, 406--413.
##' @examples
##' nlrob.tau( density ~ Asym/(1 + exp(( xmid -log(conc) )/scal )),
##'            data = DNase[DNase$Run == 1, ],
##'            pnames = c("Asym", "xmid", "scal"),
##'            lower = 0, upper = 3,  trace = TRUE )
##'
##' @author Eduardo L. T. Conceicao
nlrob.tau <- function(formula, data, pnames, lower, upper,
                      tolerance = 1e-5, f0 = NULL,
                      psi = c("bisquare", "optimal"),
                      tuning.chi.scale = NULL, tuning.chi.tau = NULL, ...)
{
    psi <- match.arg(psi)
    if (is.null(tuning.chi.scale))
        tuning.chi.scale <- switch(psi, bisquare = list(b = 0.20, cc = 1.55),
                                   optimal = list(b = 0.5, cc = 0.405))
    if (is.null(tuning.chi.tau))
        tuning.chi.tau <- switch(psi, bisquare = list(b = 0.46, cc = 6.04),
                                 optimal = list(b = 0.128, cc = 1.060))

    b1 <- tuning.chi.scale$b
    c1 <- tuning.chi.scale$cc
    b2 <- tuning.chi.tau$b
    c2 <- tuning.chi.tau$cc

    ## Preliminary psi-specific checks / computations:
    switch(psi, "bisquare" = {
	if(!require("PolynomF"))
	    stop(gettextf(
		"You must install package 'PolynomF' before you can use %s",
		"nlrob.tau(*, psi=bisquare)"), domain=NA)

	b1 <- b1/MrhoInf(c1, psi)
	b2 <- b2/MrhoInf(c2, psi)
    })

    rho1 <- function(t) Mchi(t, c1, psi)
    rho2 <- function(t) Mchi(t, c2, psi)

    rho.inv <- switch(psi, bisquare = function(y) {
        x <- polynom()
        p <- 3*x^2 - 3*x^4 + x^6
        t <- solve(p, y)
        t <- Re( t[Im(t) == 0] )
        unique( t[t >= 0] ) * c1
    }, optimal = function(y) {
        ## Salibian-Barrera, Matias, Willems, Gert, and Zamar, Ruben (2008).
        ## The fast-tau estimator for regression.
        ## Journal of Computational and Graphical Statistics 17, 659-682.
        sqrt(y/1.38) * c1 * 3
    })

    M_scale <- function(sigma, u) sum( rho1(u/sigma) )/nobs - b1
    tau_scale2 <- function(u, sigma) sigma^2 * 1/b2*sum( rho2(u/sigma) )/nobs

    objective <- function(par) {
        names(par) <- pnames
        fit <- eval( formula[[3L]], c(data, par) )
        res <- y - fit
        ## Rousseeuw, Peter J., and Leroy, Annick M. (1987).
        ## Robust Regression & Outlier Detection.
        ## John Wiley & Sons, New York, p. 137.
        med_abs_res <- median(abs(res))
        sigma <- uniroot( M_scale,
                          lower = constant[1L] * med_abs_res,
                          upper = constant[2L] * med_abs_res,
                          u = res )$root
        tau_scale2(res, sigma)
    }

    formula <- as.formula(formula)
    dataName <- substitute(data)
    varNames <- all.vars(formula)
    obsNames <- rownames(data <- as.data.frame(data))
    data <- as.list(data)# to be used as such
    if (length(formula) == 2L) { ## as nls
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 0
    }

    if (!is.character(pnames) || any(is.na(match(pnames, varNames))))
        stop("'pnames' must be a character vector named with parameters in 'formula'")
    npar <- length(pnames)
    if (length(lower) == 1)
        lower <- rep(lower, npar)
    if (length(lower) != npar)
        stop(gettextf("lower must be either of length %d, or length 1", npar))
    if (length(upper) == 1)
        upper <- rep(upper, npar)
    if (length(upper) != npar)
        stop(gettextf("upper must be either of length %d, or length 1", npar))
    stopifnot(is.numeric(lower), is.numeric(upper), lower <= upper)
    y <- eval(formula[[2L]], data)
    nobs <- length(y)
    stopifnot(nobs >= npar)
    if (is.null(f0))
        f0 <- mean(scale(y, scale = FALSE)^2)
    stopifnot(is.numeric(f0), f0 > 0)
    constant <- c(
        switch(psi,
               bisquare = 1/c1,
               optimal = 1/c1 * 1/3),
        if (nobs %% 2) 2/rho.inv(2/(nobs+2)) else 1/rho.inv(1/(nobs+1)))
    optRes <-
        JDEoptim(lower, upper, objective, tol = tolerance, fnscale = f0, ...)
    coef <- optRes$par
    names(coef) <- pnames
    iter <- optRes$iter
    status <-
        if (optRes$convergence == 0)
            "converged"
        else paste("failed to converge in", iter, "steps")
    fit <- eval( formula[[3L]], c(data, coef) )
    names(fit) <- obsNames
    structure(list(call = match.call(), formula=formula, nobs=nobs,
                   coefficients = coef,
                   fitted.values = fit,
                   residuals = y - fit,
                   crit = optRes$value,
                   Scale = sqrt(optRes$value),
                   status = status, iter = iter, data = dataName),
              class = "nlrob")
} ## nlrob.tau


##' CM-estimator for (constrained) nonlinear regression
##'
##' Copyright 2013, Eduardo L. T. Conceicao
##' Available under the GPL (>= 2)
##' @title CM-estimator of Nonlinear Regression
##' @param formula
##' @param data
##' @param pnames
##' @param lower
##' @param upper
##' @param tolerance
##' @param f0
##' @param psi
##' @param tuning.chi
##' @param ... optional arguments for optimization, e.g., \code{trace = TRUE}, passed to \code{\link{JDEoptim}(.)}.
##' @return
##' @references
##'   Mendes, B.V.M., and Tyler, D.E. (1996)
##'   Constrained M-estimation for regression.
##'   In: \emph{Robust Statistics, Data Analysis and Computer Intensive Methods},
##'   Lecture Notes in Statistics 109, Springer, New York, 299--320.
##' @examples
##'  nlrob.CM( density ~ Asym/(1 + exp(( xmid -log(conc) )/scal )),
##'            data = DNase[DNase$Run == 1, ],
##'            pnames = c("Asym", "xmid", "scal", "sigma"),
##'            lower = 0, upper = c(3, 3, 3, 0.1),
##'            trace = TRUE )
##' @author Eduardo L. T. Conceicao
nlrob.CM <- function(formula, data, pnames, lower, upper,
		     tolerance = 1e-5, f0 = NULL,
		     psi = c("bisquare", "lqq", "welsh", "optimal", "hampel", "ggw"),
		     tuning.chi = NULL, ...)
{
    psi <- match.arg(psi)
    if (is.null(tuning.chi))
        tuning.chi <- switch(psi, bisquare = list(b = 0.5, cc = 1, c = 4.835),
                             stop("unable to find constants for psi function"))

    b <- tuning.chi$b
    cc <- tuning.chi$cc
    c <- tuning.chi$c

    rho <- function(t) Mchi(t, cc, psi)
    M_scale <- function(sigma, u) sum( rho(u/sigma) )/nobs - b

    if ("sigma" %in% pnames) {
        objective <- function(par) {
            names(par) <- pnames
            fit <- eval( formula[[3L]], c(as.list(data), par) )
            sigma <- par["sigma"]
            c * sum(rho( (y - fit)/sigma ))/nobs + log(sigma)
        }
        con <- function(par) {
            names(par) <- pnames
            fit <- eval( formula[[3L]], c(as.list(data), par) )
            M_scale(par["sigma"], y - fit)
        }
    } else {
        objective <- function(par) {
            names(par) <- pnames
            fit <- eval( formula[[3L]], c(as.list(data), par) )
            resid <- y - fit
            sigma <- mad(resid)
            c * sum(rho( resid/sigma ))/nobs + log(sigma)
        }
        con <- NULL
    }


    formula <- as.formula(formula)
    dataName <- substitute(data)
    varNames <- all.vars(formula)
    data <- as.data.frame(data)
    if (length(formula) == 2L) {
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 0
    }

    if (!is.character(pnames) || any(is.na(match(pnames[!(pnames == "sigma")], varNames))))
        stop("'pnames' must be a character vector named with parameters in 'formula'")
    npar <- length(pnames)
    if (length(lower) == 1)
        lower <- rep(lower, npar)
    if (length(lower) != npar)
        stop(gettextf("lower must be either of length %d, or length 1", npar))
    if (length(upper) == 1)
        upper <- rep(upper, npar)
    if (length(upper) != npar)
        stop(gettextf("upper must be either of length %d, or length 1", npar))
    stopifnot(is.numeric(lower), is.numeric(upper), lower <= upper)
    if ("sigma" %in% pnames) {
        if ("sigma" %in% varNames || "sigma" %in% names(data))
            stop("Do not use 'sigma' as a variable name or as a parameter name in 'formula'")
        stopifnot(lower[pnames == "sigma"] >= 0)
    }
    y <- eval(formula[[2L]], as.list(data))
    nobs <- length(y)
    stopifnot(nobs >= npar)
    if (is.null(f0))
        f0 <- mean(scale(y, scale = FALSE)^2)
    stopifnot(is.numeric(f0), f0 > 0)

    optRes <- JDEoptim(lower, upper, objective, con,
                       tol = tolerance, fnscale = f0, ...)
    coef <- optRes$par
    names(coef) <- pnames
    iter <- optRes$iter
    status <- if (optRes$convergence == 0)
        "converged"
    else paste("failed to converge in", iter, "steps")
    fit <- eval( formula[[3L]], c(as.list(data), coef) )
    names(fit) <- rownames(data)

    structure(list(call = match.call(), formula=formula, nobs=nobs,
                   coefficients = coef,
                   fitted.values = fit,
                   residuals = y - fit,
                   crit = optRes$value,
                   status = status, iter = iter, data = dataName),
              class = "nlrob")
} ## nlrob.CM


##' Compute a  Maximum Trimmed Likelihood (MTL) estimator for nonlinear
##' (constrained) regression.
##'
##' Copyright 2013, Eduardo L. T. Conceicao
##' Available under the GPL (>= 2)
##' @title Maximum Trimmed Likelihood (MTL) Estimator for Nonlinear Regression
##' @param formula
##' @param data
##' @param pnames
##' @param lower
##' @param upper
##' @param cutoff
##' @param tolerance
##' @param f0
##' @param ... optional arguments for optimization, e.g., \code{trace = TRUE}, passed to \code{\link{JDEoptim}(.)}.
##' @return
##' @references
##' Hadi, Ali S., and Luceno, Alberto (1997).
##'    Maximum trimmed likelihood estimators: a unified approach,
##'    examples, and algorithms.
##'    Computational Statistics & Data Analysis 25, 251-272.
##'
##' Gervini, Daniel, and Yohai, Victor J. (2002).
##'    A class of robust and fully efficient regression estimators.
##'    The Annals of Statistics 30, 583-616.
##' @source
##' Maronna, Ricardo A., Martin, R. Douglas, and Yohai, Victor J. (2006).
##' \emph{Robust Statistics: Theory and Methods}
##' Wiley, Chichester, p. 133.
##' @examples
##' nlrob.mtl( density ~ Asym/(1 + exp(( xmid -log(conc) )/scal )),
##'            data = DNase[DNase$Run == 1, ],
##'            pnames = c("Asym", "xmid", "scal", "sigma"),
##'            lower = 0, upper = c(3, 3, 3, 0.1),
##'            trace = TRUE )
##' @author Eduardo L. T. Conceicao
nlrob.mtl <- function(formula, data, pnames, lower, upper,
                      cutoff = 2.5,
                      tolerance = 1e-5, f0 = NULL, ...)
{
    trim <- function(t) {
        t <- sort.int(t)
        i <- which(t >= cutoff)
        if (length(i) > 0) {
            h <- floor(min( (i - 1)/(2*pnorm(t[i]) - 1) ))
            h <- max(h, hlow)
        } else h <- nobs
        list(h = h, t = t)
    }

    constant <- log(2*pi)
    if ("sigma" %in% pnames) {
        objective <- function(par) {
            names(par) <- pnames
            fit <- eval( formula[[3L]], c(as.list(data), par) )
            sigma <- par["sigma"]
            tp <- trim( abs( (y - fit)/sigma ) )
            h <- tp$h
            h*(constant + 2*log(sigma)) + sum(tp$t[1L:h]^2)
        }
    } else {
        objective <- function(par) {
            names(par) <- pnames
            fit <- eval( formula[[3L]], c(as.list(data), par) )
            resid <- y - fit
            sigma <- mad(resid)
            tp <- trim( abs(resid/sigma) )
            h <- tp$h
            h*(constant + 2*log(sigma)) + sum(tp$t[1L:h]^2)
        }
    }


    formula <- as.formula(formula)
    dataName <- substitute(data)
    varNames <- all.vars(formula)
    data <- as.data.frame(data)
    if (length(formula) == 2L) {
        formula[[3L]] <- formula[[2L]]
        formula[[2L]] <- 0
    }

    if (!is.character(pnames) || any(is.na(match(pnames[!(pnames == "sigma")], varNames))))
        stop("'pnames' must be a character vector named with parameters in 'formula'")
    npar <- length(pnames)
    if (length(lower) == 1)
        lower <- rep(lower, npar)
    if (length(lower) != npar)
        stop(gettextf("lower must be either of length %d, or length 1", npar))
    if (length(upper) == 1)
        upper <- rep(upper, npar)
    if (length(upper) != npar)
        stop(gettextf("upper must be either of length %d, or length 1", npar))
    stopifnot(is.numeric(lower), is.numeric(upper), lower <= upper)
    if ("sigma" %in% pnames) {
        if ("sigma" %in% varNames || "sigma" %in% names(data))
            stop("Do not use 'sigma' as a variable name or as a parameter name in 'formula'")
        stopifnot(lower[pnames == "sigma"] >= 0)
    }
    y <- eval(formula[[2L]], as.list(data))
    nobs <- length(y)
    stopifnot(nobs >= npar)
    if (is.null(f0))
        f0 <- sum(scale(y, scale = FALSE)^2)
    stopifnot(is.numeric(f0), f0 > 0)
    hlow <- (nobs + npar + 1)%/%2

    optRes <-
        JDEoptim(lower, upper, objective, tol = tolerance, fnscale = f0, ...)
    coef <- optRes$par
    names(coef) <- pnames
    crit <- optRes$value
    iter <- optRes$iter
    status <- if (optRes$convergence == 0)
        "converged"
    else paste("failed to converge in", iter, "steps")
    fit <- eval( formula[[3L]], c(as.list(data), coef) )
    names(fit) <- rownames(data)
    resid <- y - fit
    quan <-
        trim( resid/(if ("sigma" %in% pnames) coef["sigma"] else mad(resid)) )$h

    structure(list(call = match.call(), formula=formula, nobs=nobs,
                   coefficients = coef,
                   fitted.values = fit,
                   residuals = resid,
                   crit = crit,
                   quan = quan,
                   status = status, iter = iter, data = dataName),
              class = "nlrob")
} ## nlrob.mtl
