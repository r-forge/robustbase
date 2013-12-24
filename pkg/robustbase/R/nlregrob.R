### TODOs / Questions (Martin  <-->  Eduardo):
### -------------------------------------------

## (2) argument names:
## 5 matches for "@param .*tuning" in buffer  nlregrob.R
##      27:##' @param tuning.chi.scale \
##      28:##' @param tuning.psi.M     /  nlrob.MM
##     217:##' @param tuning.chi.scale \
##     218:##' @param tuning.chi.tau   /  nlrob.tau
##     364:##' @param tuning.chi          nlrob.CM


nlrob.MM <-
    function(formula, data, pnames, lower, upper, tol = 1e-6,
	     psi = c("bisquare", "lqq", "optimal", "hampel"),
	     init = c("S", "lts"),
             ctrl = nlrob.control("MM", psi=psi, init=init, fnscale=NULL,
                 tuning.chi.scale = .psi.conv.cc(psi, .Mchi.tuning.defaults[[psi]]),
                 tuning.psi.M     = .psi.conv.cc(psi, .Mpsi.tuning.defaults[[psi]]),
                 optim.control = list(), jOptimArgs = list(...)),
             ...)
{
    if(missing(ctrl)) {
        init <- match.arg(init)
        psi <- match.arg(psi)
        force(ctrl) #
    } else {
        init <- ctrl$ init
        psi  <- ctrl$ psi
    }
    c1 <- ctrl$tuning.chi.scale
    c2 <- ctrl$tuning.psi.M

    ## Preliminary psi-specific checks / computations:
    switch(psi,
           "lqq" = { # lqqMax = rho(Inf), used in rho.inv() *and* 'constant':
               c12 <- c1[1]+c1[2]
               lqqMax <- (c1[1]*c1[3] - 2*c12)/(1-c1[3]) + c12})

    rho1 <- function(t) Mchi(t, c1, psi)
    rho2 <- function(t) Mchi(t, c2, psi)

    rho.inv <- switch(psi, "bisquare" = function(y) {
        ## Find  x := u^2  which solves cubic eq. 3*x - 3*x^2 + x^3 = y
        ##      <==> (x-1)^3 + 1 = y  <==> (1-x)^3 = 1-y  <==> x = 1 - (1-y)^(1/3)
        ##      (where we assume 0 <= y <= 1, i.e, y-1 < 0)
        c1 * sqrt(1 - (1 - y)^(1/3))
    }, "lqq" = function(y) {
        uniroot( function(x) rho1(x) - y, lower = 0, upper = lqqMax )$root
    }, "optimal" = function(y) {
        ## Salibian-Barrera, Matias, Willems, Gert, and Zamar, Ruben (2008).
        ## The fast-tau estimator for regression.
        ## Journal of Computational and Graphical Statistics 17, 659-682.
        sqrt(y/1.38) * c1 * 3
    }, "hampel" = function(y) {
        C <- MrhoInf(c1, psi)
        a <- c1[1]; b <- c1[2]; r <- c1[3]
        if (y <= a/C)
            sqrt(2*C*y)
        else if (y <= (2*b - a)/C)
            0.5*a + C/a*y
        else r + sqrt( r^2 - ( (r - b)*(2*C/a*y + (b - a)) - b*r ) )
    }, stop(gettextf("Psi function '%s' not supported yet", psi)))

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
	       }, stop(gettextf("Initialization 'init = \"%s\"' not supported (yet)",
				init)))

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
    if (is.null(fnscale <- ctrl$ fnscale))
        fnscale <- sum((y - mean(y))^2)
    ctrl$fnscale <- NULL # remove it there
    stopifnot(is.numeric(fnscale), fnscale > 0)
    constant <- c(
        switch(psi, bisquare = 1/c1,
               lqq = 1/lqqMax,
               optimal = 1/c1 * 1/3,
               hampel = 1/c1[3]),
        if(nobs %% 2) 2/rho.inv(2/(nobs+2)) else 1/rho.inv(1/(nobs+1)))
    switch(init, lts = h <- (nobs + npar + 1)%/%2)

    initial <- do.call(JDEoptim, c(list(lower, upper, objective.initial,
					tol=tol, fnscale=fnscale), ctrl$jOptimArgs))
    names(initial$par) <- pnames
    res <- y - eval( formula[[3L]], c(data, initial$par) )

    med_abs_res <- median(abs(res))
    sigma <- uniroot( M_scale,
                      lower = constant[1L] * med_abs_res,
                      upper = constant[2L] * med_abs_res,
                      u = res )$root

    M <- optim(initial$par, objective.M, sigma = sigma,
               method = "L-BFGS-B", lower = lower, upper = upper,
               control = c(list(fnscale = initial$value, parscale = initial$par),
                           ctrl$optim.control) )
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
                   status = status, counts = M$counts, data = dataName, ctrl=ctrl),
              class = "nlrob")
} ## nlrob.MM


nlrob.tau <- function(formula, data, pnames, lower, upper, tol = 1e-6,
		      psi = c("bisquare", "optimal"),
		      ctrl = nlrob.control("tau", psi=psi, fnscale=NULL,
			  tuning.chi.scale = NULL, tuning.chi.tau = NULL,
			  jOptimArgs = list(...)),
		      ...)
{
    if(missing(ctrl)) {
	psi <- match.arg(psi)
	force(ctrl) #
    } else {
	psi  <- ctrl$ psi
    }
    if(is.null(.chi.s <- ctrl$tuning.chi.scale))
	.chi.s <- switch(psi, bisquare = list(b = 0.20, cc = 1.55),
			 optimal = list(b = 0.5, cc = 0.405))
    if(is.null(.chi.t <- ctrl$tuning.chi.tau))
	.chi.t <- switch(psi, bisquare = list(b = 0.46, cc = 6.04),
			 optimal = list(b = 0.128, cc = 1.060))

    b1 <- .chi.s$b
    c1 <- .chi.s$cc
    b2 <- .chi.t$b
    c2 <- .chi.t$cc

    ## Preliminary psi-specific checks / computations:
    switch(psi, "bisquare" = {
	b1 <- b1/MrhoInf(c1, psi)
	b2 <- b2/MrhoInf(c2, psi)
    })

    rho1 <- function(t) Mchi(t, c1, psi)
    rho2 <- function(t) Mchi(t, c2, psi)

    rho.inv <- switch(psi, "bisquare" = function(y) {
	c1 * sqrt(1 - (1 - y)^(1/3))
    }, "optimal" = function(y) {
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
    if (is.null(fnscale <- ctrl$ fnscale))
        fnscale <- mean((y - mean(y))^2)
    ctrl$fnscale <- NULL # remove it there
    stopifnot(is.numeric(fnscale), fnscale > 0)
    constant <- c(
        switch(psi,
               bisquare = 1/c1,
               optimal = 1/c1 * 1/3),
        if (nobs %% 2) 2/rho.inv(2/(nobs+2)) else 1/rho.inv(1/(nobs+1)))

    optRes <- do.call(JDEoptim, c(list(lower, upper, objective, tol=tol, fnscale=fnscale),
				  ctrl$jOptimArgs))
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
                   status = status, iter = iter, data = dataName, ctrl=ctrl),
              class = "nlrob")
} ## nlrob.tau


nlrob.CM <- function(formula, data, pnames, lower, upper, tol = 1e-6,
		     psi = c("bisquare", "lqq", "welsh", "optimal", "hampel", "ggw"),
                     ctrl = nlrob.control("CM", psi=psi, fnscale=NULL,
                         tuning.chi = NULL, jOptimArgs = list(...)),
                     ...)
{
    if(missing(ctrl)) {
        psi <- match.arg(psi)
        force(ctrl) #
    } else {
        psi  <- ctrl$ psi
    }
    if (is.null(t.chi <- ctrl$tuning.chi))
	t.chi <- switch(psi, bisquare = list(b = 0.5, cc = 1, c = 4.835),
			stop("unable to find constants for psi function"))
    b  <- t.chi$b
    cc <- t.chi$cc
    c  <- t.chi$c

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
    if (is.null(fnscale <- ctrl$ fnscale))
        fnscale <- mean((y - mean(y))^2)
    ctrl$fnscale <- NULL # remove it there
    stopifnot(is.numeric(fnscale), fnscale > 0)

    optRes <- do.call(JDEoptim, c(list(lower, upper, objective, constr=con,
				       tol=tol, fnscale=fnscale),
				  ctrl$jOptimArgs))
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
                   status = status, iter = iter, data = dataName, ctrl=ctrl),
              class = "nlrob")
} ## nlrob.CM


nlrob.mtl <- function(formula, data, pnames, lower, upper, tol = 1e-6,
                      ctrl = nlrob.control("mtl", cutoff = 2.5, jOptimArgs = list(...)),
                      ...)
{
    cutoff <- ctrl[["cutoff"]]
    trim <- function(t) {
        t <- sort.int(t)
        i <- which(t >= cutoff)
        h <- if (length(i) > 0)
            max(hlow, floor(min( (i - 1)/(2*pnorm(t[i]) - 1) ))) else nobs
        list(h = h, t = t)
    }

    constant <- log(2*pi)
    if ("sigma" %in% pnames) {
        objective <- function(par) {
            names(par) <- pnames
            fit <- eval( formula[[3L]], c(as.list(data), par) )
            sigma <- par[["sigma"]]
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
    if (is.null(fnscale <- ctrl$ fnscale))
        fnscale <- sum((y - mean(y))^2)
    ctrl$fnscale <- NULL # remove it there
    stopifnot(is.numeric(fnscale), fnscale > 0)
    hlow <- (nobs + npar + 1)%/%2

    optRes <- do.call(JDEoptim, c(list(lower, upper, objective, tol=tol, fnscale=fnscale),
				  ctrl$jOptimArgs))
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

    structure(list(call = match.call(), formula=formula, nobs=nobs, method="mtl",
                   coefficients = coef,
                   fitted.values = fit,
                   residuals = resid,
                   crit = crit,
                   quan = quan,
                   status = status, iter = iter, data = dataName, ctrl = ctrl),
              class = "nlrob")
} ## nlrob.mtl

nlrob.control <- function(method,
                          psi = c("bisquare", "lqq", "welsh", "optimal", "hampel", "ggw"),
                          init = c("S", "lts"),
                          jOptimArgs  = list(),
                          ...)
{
  psi <- match.arg(psi)
  init <- match.arg(init)
  dots <- list(...)
  argNms <- names(dots)
  ##' argument or default -> return list of length 1
  a. <- function(nm,def) {
    L <- list( if(nm %in% argNms) dots[[nm]] else def )
    names(L) <- nm
    L
  }
  switch(method,
         "M" = {
             list(method = method) # not yet used
         },
         "MM" = {
             c(list(method = method, init = init, psi = psi),
               a.("fnscale", NULL),
               a.("tuning.chi.scale", .psi.conv.cc(psi, .Mchi.tuning.defaults[[psi]])),
               a.("tuning.psi.M",     .psi.conv.cc(psi, .Mpsi.tuning.defaults[[psi]])),
               a.("optim.control", list()),
               list(jOptimArgs = jOptimArgs))
         },
         "tau" = {
             c(list(method = method, psi = psi),
               a.("fnscale", NULL),
               a.("tuning.chi.scale", NULL),
               a.("tuning.chi.tau", NULL),
               list(jOptimArgs = jOptimArgs))
         },
         "CM" = {
             c(list(method = method, psi = psi),
               a.("fnscale", NULL),
               a.("tuning.chi", NULL),
               list(jOptimArgs = jOptimArgs))
         },
         "mtl" = {
             c(list(method = method),
               a.("fnscale", NULL),
               a.("cutoff", 2.5),
               list(jOptimArgs = jOptimArgs))
         },
         stop("Method ", method, "not correctly supported yet"))
}
