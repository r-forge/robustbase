##' Nonlinear constrained minimization via Differential Evolution (jDE) [1].
##' It uses the DE/rand/1/either-or mutation strategy [2].
##' It also utilizes only one population [3] against two sets as in the
##' original DE algorithm [4].
##'
##' Copyright 2013, Eduardo L. T. Conceicao
##' Available under the GPL (>= 2)
##' @title Nonlinear Constrained Minimization via Differential Evolution
##' @param lower
##' @param upper
##' @param fn
##' @param constr
##' @param meq
##' @param eps
##' @param NP
##' @param Fl
##' @param Fu
##' @param tau1
##' @param tau2
##' @param tau3
##' @param jitter_factor
##' @param tol
##' @param maxiter
##' @param fnscale
##' @param FUN
##' @param add_to_init_pop
##' @param trace
##' @param triter
##' @param verbose
##' @param ... optional arguments passed to \code{fn()} \emph{and}
##' \code{constr()} if that is not \code{NULL}
##' @return
##' @seealso Function \code{\link{DEoptim}{DEoptim}} in the \pkg{DEoptim} package
##'  has many more options than \code{JDEoptim()}, but does not allow constraints
##' in the same flexible manner.
##' @references
##' [1] Brest, J., Greiner, S., Boskovic, B., Mernik, M. and Zumer, V. (2006).
##'     Self-adapting control parameters in differential evolution: a comparative
##'     study on numerical benchmark problems.
##'     IEEE Trans. Evol. Comput. 10, 646-657.
##'
##' [2] Price, Kenneth V., Storn, Rainer M., and Lampinen, Jouni A. (2005).
##'     Differential Evolution: a practical approach to global optimization.
##'     Springer, Berlin, pp. 117-118.
##'
##' [3] Babu, B. V., and Angira, Rakesh (2006).
##'     Modified differential evolution (MDE) for optimization of non-linear
##'     chemical processes.
##'     Computers and Chemical Engineering 30, 989-1002.
##'
##' [4] Storn, Rainer, and Price, Kenneth (1997).
##'     Differential evolution - a simple and efficient heuristic for global
##'     optimization over continuous spaces.
##'     Journal of Global Optimization 11, 341-359.
##' @examples
##' JDEoptim(c(-100, -100), c(100, 100), sf1, tol = 1e-7, trace = TRUE)
##' JDEoptim(rep(-500, 10), rep(500, 10), swf, tol = 1e-7, trace = TRUE)
##'
##' JDEoptim(c(1e-5, 1e-5), c(16, 16),
##'          RND$obj, RND$con, tol = 1e-7, trace = TRUE)
##' JDEoptim(c(100, 1000, 1000, 10, 10), c(10000, 10000, 10000, 1000, 1000),
##'          HEND$obj, HEND$con, tol = 1e-7, trace = TRUE)
##' JDEoptim(c(1500, 1, 3000, 85, 90, 3, 145), c(2000, 120, 3500, 93, 95, 12, 162),
##'          alkylation$obj, alkylation$con, tol = 1e-7, trace = TRUE)
##' @author Eduardo L. T. Conceicao
JDEoptim <-
    function(lower, upper, fn, constr = NULL, meq = 0, eps = 1e-5,
             NP = 10*d, Fl = 0.1, Fu = 1, tau1 = 0.1, tau2 = 0.1, tau3 = 0.1,
             jitter_factor = 0.001,
             tol = 1e-15, maxiter = 200*d, fnscale = 1,
             FUN = c("median", "max"),
             add_to_init_pop = NULL, trace = FALSE, triter = 1,
             details = FALSE, ...)
{
    handle.bounds <- function(x, u) {
    # Check feasibility of bounds and enforce parameters limits
    # by a deterministic variant of bounce-back resetting
    # Price, KV, Storn, RM, and Lampinen, JA (2005) Differential Evolution.
    # Springer, p 206
        bad <- x > upper
        x[bad] <- 0.5*(upper[bad] + u[bad])
        bad <- x < lower
        x[bad] <- 0.5*(lower[bad] + u[bad])
        x
    }

    performReproduction <- function() {
        ignore <- runif(d) > CRtrial
        if (all(ignore))                  # ensure that trial gets at least
            ignore[sample(d, 1)] <- FALSE # one mutant parameter
        # Source for trial is the base vector plus weighted differential
        trial <- if (runif(1) <= pF[i])
            X.base + Ftrial*(X.r1 - X.r2)
        else X.base + 0.5*(Ftrial + 1)*(X.r1 + X.r2 - 2*X.base)
        # or trial parameter comes from target vector X.i itself.
        trial[ignore] <- X.i[ignore]
        trial
    }

    child <- if (is.null(constr)) {
        expression({
            ftrial <- fn1(trial) # Evaluate trial with your function
            if (ftrial <= fpop[i]) {
                pop[, i] <- trial
                fpop[i] <- ftrial
                F[, i] <- Ftrial
                CR[i] <- CRtrial
            }
        })
    } else {
        # Zhang, Haibo, and Rangaiah, G. P. (2012).
        # An efficient constraint handling method with integrated differential
        # evolution for numerical and engineering optimization.
        # Computers and Chemical Engineering 37, 74-88.
        expression({
            htrial <- constr1(trial)
            TAVtrial <- sum( pmax(htrial, 0) )
            if (TAVtrial > mu) {
                if (TAVtrial <= TAVpop[i]) {
                    pop[, i] <- trial
                    hpop[, i] <- htrial
                    F[, i] <- Ftrial
                    CR[i] <- CRtrial
                    TAVpop[i] <- TAVtrial
                }
            } else if (TAVpop[i] > mu) {
                pop[, i] <- trial
                fpop[i] <- fn1(trial)
                hpop[, i] <- htrial
                F[, i] <- Ftrial
                CR[i] <- CRtrial
                TAVpop[i] <- TAVtrial
                FF <- sum(TAVpop <= mu)/NP
                mu <- mu*(1 - FF/NP)
            } else {
                ftrial <- fn1(trial) # Evaluate trial with your function
                if (ftrial <= fpop[i]) {
                    pop[, i] <- trial
                    fpop[i] <- ftrial
                    hpop[, i] <- htrial
                    F[, i] <- Ftrial
                    CR[i] <- CRtrial
                    TAVpop[i] <- TAVtrial
                    FF <- sum(TAVpop <= mu)/NP
                    mu <- mu*(1 - FF/NP)
                }
            }
        })
    }

    if (!is.null(constr)) {
        which.best <- function(x) {
            ind <- TAVpop <= mu
            if (all(ind))
                which.min(x)
            else if (any(ind))
                which(ind)[which.min(x[ind])]
            else which.min(TAVpop)
        }
    } else which.best <- function(x) which.min(x)


    # Check input parameters
    FUN <- match.arg(FUN)
    d <- length(lower)
    if (length(upper) != d)
        stop("'lower' must have same length as 'upper'")
    stopifnot(is.numeric(lower), is.numeric(upper), lower <= upper)
    stopifnot(length(fnscale) == 1, is.numeric(fnscale), fnscale > 0)
    stopifnot(is.function(fn))
    if (!is.null(constr)) {
        stopifnot(is.function(constr))
        stopifnot(length(meq) == 1, meq == as.integer(meq), meq >= 0)
        if (length(eps) == 1)
            eps <- rep(eps, meq)
        if (length(eps) != meq)
            stop("slack must be either of length meq, or length 1")
    }
    stopifnot(length(Fl) == 1, is.numeric(Fl),
              length(Fu) == 1, is.numeric(Fu),
              Fl <= Fu)
    stopifnot(length(tau1) == 1, is.numeric(tau1), tau1 >= 0, tau1 <= 1)
    stopifnot(length(tau2) == 1, is.numeric(tau2), tau2 >= 0, tau2 <= 1)
    stopifnot(length(tau3) == 1, is.numeric(tau3), tau3 >= 0, tau3 <= 1)
    if (!is.null(jitter_factor))
        stopifnot(length(jitter_factor) == 1, is.numeric(jitter_factor))
    stopifnot(length(tol) == 1, is.numeric(tol))
    stopifnot(length(maxiter) == 1, maxiter == as.integer(maxiter))
    stopifnot(length(triter) == 1, triter == as.integer(triter))
    if (!is.null(add_to_init_pop)) {
        stopifnot(NROW(add_to_init_pop) == d)
        stopifnot(is.numeric(add_to_init_pop),
                  add_to_init_pop >= lower,
                  add_to_init_pop <= upper)
    }

    # Set default values
    defaultopt.jitter <- if (is.null(jitter_factor)) FALSE else TRUE

    # Initialization:
    fn1 <- function(par) fn(par, ...)

    if (!is.null(constr)) {
        if (meq > 0) {
            equalIndex <- 1:meq
            constr1 <- function(par) {
                h <- constr(par, ...)
                h[equalIndex] <- abs(h[equalIndex]) - eps
                h
            }
        } else constr1 <- function(par) constr(par, ...)
    }

    # Zielinski, Karin, and Laur, Rainer (2008).
    # Stopping criteria for differential evolution in
    # constrained single-objective optimization.
    # In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
    # SCI 143, Springer-Verlag, pp 111-138
    conv <- expression(( do.call(FUN, list(fpop)) - fpop[x.best.ind] )/fnscale)
    pop <- matrix(runif(NP*d, lower, upper), nrow = d)
    if (!is.null(add_to_init_pop)) {
        pop <- unname(cbind(pop, add_to_init_pop))
        NP <- ncol(pop)
    }
    stopifnot(length(NP) == 1, NP == as.integer(NP), NP >= 4)
    F <- if (defaultopt.jitter)
        (1 + jitter_factor*runif(d, -0.5, 0.5)) %o% runif(NP, Fl, Fu)
    else matrix(runif(NP, Fl, Fu), nrow = 1)
    CR <- runif(NP)
    pF <- runif(NP)
    fpop <- apply(pop, 2, fn1)
    if (!is.null(constr)) {
        hpop <- apply(pop, 2, constr1)
        if ( any(is.na(hpop)) )
            stop("value of meq is invalid")
        if (is.vector(hpop)) dim(hpop) <- c(1, length(hpop))
        TAVpop <- apply( hpop, 2, function(x) sum(pmax(x, 0)) )
        mu <- median(TAVpop)
    }

    popIndex <- 1:NP
    x.best.ind <- which.best(fpop)
    converge <- eval(conv)
    rule <- if (!is.null(constr))
        expression(converge >= tol || any(hpop[, x.best.ind] > 0))
    else expression(converge >= tol)
    convergence <- 0
    iteration <- 0

    while (eval(rule)) { # generation loop
        iteration <- iteration + 1
        if (iteration > maxiter) {
            warning("maximum number of iterations reached without convergence")
            convergence <- 1
            break
        }

        for (i in popIndex) { # Start loop through population
            # Fi update
            # Combine jitter with dither
            # Storn, Rainer (2008).
            # Differential evolution research - trends and open questions.
            # In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
            # SCI 143, Springer-Verlag, pp 11-12
            Ftrial <- if (runif(1) <= tau1) {
                if (defaultopt.jitter)
                    runif(1, Fl, Fu) * (1 + jitter_factor*runif(d, -0.5, 0.5))
                else runif(1, Fl, Fu)
            } else F[, i]
            # CRi update
            CRtrial <- if (runif(1) <= tau2)
                runif(1)
            else CR[i]
            # pFi update
            if (runif(1) <= tau3)
                pF[i] <- runif(1)

            # DE/rand/1/either-or/bin
            X.i <- pop[, i]
            # Randomly pick 3 vectors all diferent from target vector
            r <- sample(popIndex[-i], 3)
            X.base <- pop[, r[1L]]
            X.r1 <- pop[, r[2L]]
            X.r2 <- pop[, r[3L]]

            trial <- handle.bounds(performReproduction(), X.base)

            eval(child)

            x.best.ind <- which.best(fpop)
        }

        converge <- eval(conv)
        if (trace && (iteration %% triter == 0))
            cat(iteration, ":", "<", converge, ">", "(", fpop[x.best.ind], ")",
                pop[, x.best.ind],
                if (!is.null(constr))
                    paste("{", which(hpop[, x.best.ind] > 0), "}"),
                fill = TRUE)
    }

    res <- list(par = pop[, x.best.ind],
                value = fpop[x.best.ind],
                iter = iteration - 1,
                convergence = convergence)
    if (details) {
        res$poppar <- pop
        res$popcost <- fpop
    }
    res
}

## rather the whole package via 'ByteCompile: yes' in ../DESCRIPTION ?
library(compiler)
JDEoptim <- cmpfun(JDEoptim)
