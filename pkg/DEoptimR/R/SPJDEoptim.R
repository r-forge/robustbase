SPJDEoptim <- function(
    lower, upper, fn, constr = NULL, meq = 0, eps = 1e-5,
    sequential_eval = c("none", "constraints", "objective", "both"),
    ..., further_args = list(),
    NP = 10*length(lower), Fl = 0.1, Fu = 1, CRl = 0, CRu = 1.1,
    tau_F = 0.1, tau_CR = 0.1, tau_pF = 0.1,
    jitter_factor = 0.001,
    tol = 1e-15, maxiter = 2000*length(lower), fnscale = 1,
    compare_to = c("median", "max"),
    add_to_init_pop = NULL,
    trace = FALSE, triter = 1,
    details = FALSE) {

#   Copyright 2026, Eduardo L. T. Conceicao
#   Available under the GPL (>= 2)

    handle_bounds <- function(x, u) {
        # Check feasibility of bounds and enforce parameters limits
        # by a deterministic variant of bounce-back resetting
        # (also known as midpoint target/base)
        # Price, KV, Storn, RM, and Lampinen, JA (2005)
        # Differential Evolution: A Practical Approach to Global Optimization.
        # Springer, p 206
        bad <- x > upper
        x[bad] <- 0.5*(upper[bad] + u[bad])
        bad <- x < lower
        x[bad] <- 0.5*(lower[bad] + u[bad])
        x
    }

    perform_reproduction <- function() { # Mutate/recombine
        ignore <- runif(d) > CRtrial[i]
        if (all(ignore))                  # ensure that trial gets at least
            ignore[sample(d, 1)] <- FALSE # one mutant parameter

        # Source for trial is the base vector plus weighted differential
        trial <- if (runif(1) <= pFtrial[i])
            X_base + Ftrial[, i]*(X_r1 - X_r2)
        else X_base + 0.5*(Ftrial[, i] + 1)*(X_r1 + X_r2 - 2*X_base)

        # or trial parameter comes from target vector X_i itself.
        trial[ignore] <- X_i[ignore]
        trial
    }

    mirai_enabled <- requireNamespace("mirai", quietly = TRUE)

    # Check input parameters
    sequential_eval <- match.arg(sequential_eval)
    compare_to <- match.arg(compare_to)
    stopifnot(
        length(upper) == length(lower),
        length(lower) > 0, is.numeric(lower), is.finite(lower),
        length(upper) > 0, is.numeric(upper), is.finite(upper),
        lower <= upper,
        is.function(fn), length(formals(fn)) > 0
    )
    if (!is.null(constr))
        stopifnot(
            is.function(constr),
            length(formals(constr)) == length(formals(fn)),
            identical( formalArgs(constr)[-1L], formalArgs(fn)[-1L] ),
            length(meq) == 1, meq == as.integer(meq), meq >= 0,
            is.numeric(eps), is.finite(eps), eps > 0,
            length(eps) == 1 || length(eps) == meq
        )
    stopifnot(is.list(further_args))
    if (length(further_args) > 0) {
        stopifnot(
            length(names(further_args)) > 0,
            all( nzchar(names(further_args)) )
        )
        if ( !("..." %in% formalArgs(fn)) )
            stopifnot(
                all( names(further_args) %in% formalArgs(fn)[-1L] ),
                all( formalArgs(fn)[-1L][sapply(formals(fn)[-1L], is.name)] %in%
                     names(further_args) )
            )
    }
    stopifnot(
        length(NP) == 1, NP == as.integer(NP), NP >= 0,
        length(Fl) == 1, is.numeric(Fl), is.finite(Fl),
        length(Fu) == 1, is.numeric(Fu), is.finite(Fu),
        Fl <= Fu,
        length(CRl) == 1, is.numeric(CRl), CRl >= 0, CRl <= 1,
        length(CRu) == 1, is.numeric(CRu), CRu >= 0,
        CRl <= CRu
    )
    stopifnot(
        length(tau_F) == 1, is.numeric(tau_F), tau_F >= 0, tau_F <= 1,
        length(tau_CR) == 1, is.numeric(tau_CR), tau_CR >= 0, tau_CR <= 1,
        length(tau_pF) == 1, is.numeric(tau_pF), tau_pF >= 0, tau_pF <= 1
    )
    if (!is.null(jitter_factor))
        stopifnot(
            length(jitter_factor) == 1,
            is.numeric(jitter_factor),
            is.finite(jitter_factor)
        )
    stopifnot(
        length(tol) == 1, is.numeric(tol), is.finite(tol),
        length(maxiter) == 1, maxiter == as.integer(maxiter), maxiter >= 0,
        length(fnscale) == 1, is.numeric(fnscale),
        is.finite(fnscale), fnscale > 0
    )
    if (!is.null(add_to_init_pop))
        stopifnot(
            NROW(add_to_init_pop) == length(lower),
            is.numeric(add_to_init_pop),
            is.finite(add_to_init_pop),
            add_to_init_pop >= lower,
            add_to_init_pop <= upper
        )
    stopifnot(
        length(trace) == 1, is.logical(trace), !is.na(trace),
        length(triter) == 1, triter == as.integer(triter), triter >= 1,
        length(details) == 1, is.logical(details), !is.na(details)
    )

    has_NO_constraints <- is.null(constr)
    has_constraints <- !has_NO_constraints

    which_best <- if (has_constraints)
        function(x) {
            ind <- TAVpop <= mu
            if (all(ind))
                which.min(x)
            else if (any(ind))
                which(ind)[which.min(x[ind])]
            else which.min(TAVpop)
        }
    else which.min

    child <- if (has_NO_constraints) { # Evaluate/select
        expression({
            ftrial <- fn2(trial)
            i <- ftrial <= fpop
            pop[, i] <- trial[, i]
            fpop[i] <- ftrial[i]
            F[, i] <- Ftrial[, i]
            CR[i] <- CRtrial[i]
            pF[i] <- pFtrial[i]
        })
    } else if (meq > 0) { # equality constraints are present
                          # alongside the inequalities
        # Zhang, Haibo, and Rangaiah, G. P. (2012).
        # An efficient constraint handling method with integrated differential
        # evolution for numerical and engineering optimization.
        # Computers and Chemical Engineering 37, 74-88.
        expression({
            htrial <- constr2(trial)
            TAVtrial <- colSums( pmax(htrial, 0) )
            feas_trial <- TAVtrial <= mu
            feas_target <- TAVpop <= mu
            if ( any(i <- feas_trial) )
                ftrial[i] <- fn2(trial[, i])

            if ( any(i <- !feas_trial & TAVtrial <= TAVpop) ) {
                pop[, i] <- trial[, i]   # trial and target are both
                hpop[, i] <- htrial[, i] # infeasible, the one with smaller
                F[, i] <- Ftrial[, i]    # constraint violation is chosen
                CR[i] <- CRtrial[i]      # or trial vector when both are
                pF[i] <- pFtrial[i]      # solutions of equal quality
                TAVpop[i] <- TAVtrial[i]
            }

            if ( any(i <- feas_trial & !feas_target) ) {
                pop[, i] <- trial[, i]   # trial is feasible and target is not
                fpop[i] <- ftrial[i]
                hpop[, i] <- htrial[, i]
                F[, i] <- Ftrial[, i]
                CR[i] <- CRtrial[i]
                pF[i] <- pFtrial[i]
                TAVpop[i] <- TAVtrial[i]
            }

            if ( any(i <- feas_trial & feas_target & ftrial <= fpop) ) {
                pop[, i] <- trial[, i]   # between two feasible solutions,
                fpop[i] <- ftrial[i]     # the one with better objective
                hpop[, i] <- htrial[, i] # function value is chosen
                F[, i] <- Ftrial[, i]    # or trial vector when both are
                CR[i] <- CRtrial[i]      # solutions of equal quality
                pF[i] <- pFtrial[i]
                TAVpop[i] <- TAVtrial[i]
                FF <- sum(TAVpop <= mu)/NP
                mu <- mu*(1 - FF/NP)^sum(i)
            }
        })
    } else {              # only inequality constraints are present
        expression({
            htrial <- constr2(trial)
            TAVtrial <- colSums( pmax(htrial, 0) )
            feas_trial <- TAVtrial <= mu
            feas_target <- TAVpop <= mu
            if ( any(i <- feas_trial) )
                ftrial[i] <- fn2(trial[, i])

            if ( any(i <- !feas_trial & TAVtrial <= TAVpop) ) {
                pop[, i] <- trial[, i]   # trial and target both infeasible
                hpop[, i] <- htrial[, i]
                F[, i] <- Ftrial[, i]
                CR[i] <- CRtrial[i]
                pF[i] <- pFtrial[i]
                TAVpop[i] <- TAVtrial[i]
            }

            if ( any(i <- feas_trial & !feas_target) ) {
                pop[, i] <- trial[, i]   # trial is feasible and target is not
                fpop[i] <- ftrial[i]
                hpop[, i] <- htrial[, i]
                F[, i] <- Ftrial[, i]
                CR[i] <- CRtrial[i]
                pF[i] <- pFtrial[i]
                TAVpop[i] <- TAVtrial[i]
                FF <- sum(TAVpop <= mu)/NP
                mu <- mu*(1 - FF/NP)^sum(i)
            }

            if ( any(i <- feas_trial & feas_target & ftrial <= fpop) ) {
                pop[, i] <- trial[, i]   # two feasible solutions
                fpop[i] <- ftrial[i]
                hpop[, i] <- htrial[, i]
                F[, i] <- Ftrial[, i]
                CR[i] <- CRtrial[i]
                pF[i] <- pFtrial[i]
                TAVpop[i] <- TAVtrial[i]
                FF <- sum(TAVpop <= mu)/NP
                mu <- mu*(1 - FF/NP)^sum(i)
            }
        })
    }

    fn1 <- function(par) do.call(fn, c(list(par), further_args))

    fn2 <- if (sequential_eval %in% c("objective", "both"))
        function(pop) if (NCOL(pop) > 1L) apply(pop, 2, fn1) else fn1(pop)
    else {
        if (!mirai_enabled) stop("the mirai package is required")
        function(pop) {
            if (NCOL(pop) > 1L)
                mirai::mirai_map(
                    lapply(seq_len(ncol(pop)), function(j) pop[, j]), fn1, ...
                )[mirai::.stop, mirai::.flat]
            else fn1(pop)
        }
    }

    if (has_constraints) {
        constr1 <- function(par) do.call(constr, c(list(par), further_args))

        constr2 <- if (sequential_eval %in% c("constraints", "both"))
            function(pop) matrix(apply(pop, 2, constr1), ncol = NP)
        else {
            if (!mirai_enabled) stop("the mirai package is required")
            function(pop) {
                matrix(
                    mirai::mirai_map(
                        lapply(seq_len(NP), function(j) pop[, j]), constr1, ...
                    )[mirai::.stop, mirai::.flat],
                    ncol = NP
                )
            }
        }
    }

    use_jitter <- !is.null(jitter_factor)

    # Zielinski, Karin, and Laur, Rainer (2008).
    # Stopping criteria for differential evolution in
    # constrained single-objective optimization.
    # In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
    # SCI 143, Springer-Verlag, pp 111-138
    conv <- expression(
      ( do.call(compare_to, list(fpop)) - fpop[x_best_ind] )/fnscale
    )

    # Initialization
    d <- length(lower)
    pop <- matrix(runif(NP*d, lower, upper), nrow = d)
    if (!is.null(add_to_init_pop)) {
        pop <- unname(cbind(pop, add_to_init_pop))
        NP <- ncol(pop)
    }
    stopifnot(NP >= 4)
    # Combine jitter with dither
    # Storn, Rainer (2008).
    # Differential evolution research - trends and open questions.
    # In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
    # SCI 143, Springer-Verlag, pp 11-12
    F <- if (use_jitter)
        (1 + jitter_factor*runif(d, -0.5, 0.5)) %o% runif(NP, Fl, Fu)
    else matrix(runif(NP, Fl, Fu), nrow = 1)
    CR <- runif(NP, CRl, CRu)
    pF <- runif(NP)
    trial <- pop
    Ftrial <- F
    CRtrial <- CR
    pFtrial <- pF
    if (has_constraints) {
        hpop <- constr2(pop)
        stopifnot(
            is.matrix(hpop),
            !is.na(hpop), !is.nan(hpop), !is.logical(hpop),
            nrow(hpop) >= meq
        )
        if (meq > 0) {
            equal_index <- 1:meq
            constr1 <- function(par) {
                h <- do.call(constr, c(list(par), further_args))
                h[equal_index] <- abs(h[equal_index]) - eps
                h
            }
            hpop[equal_index, ] <- abs(hpop[equal_index, ]) - eps
        }
        TAVpop <- apply( hpop, 2, function(x) sum(pmax(x, 0)) )
        mu <- median(TAVpop)
        fpop <- rep_len(Inf, NP)
        if ( any(i <- TAVpop <= mu) )
            fpop[i] <- fn2(pop[, i])
        ftrial <- fpop
    } else fpop <- fn2(pop)
    stopifnot(
        is.vector(fpop),
        !is.na(fpop), !is.nan(fpop), !is.logical(fpop)
    )

    pop_index <- 1:NP
    x_best_ind <- which_best(fpop)
    converge <- eval(conv)
    rule <- if (has_constraints)
        expression(converge >= tol || any(hpop[, x_best_ind] > 0))
    else expression(converge >= tol)
    convergence <- 0
    iteration <- 0

    while (eval(rule)) { # Generation loop
        if (iteration >= maxiter) {
            warning("maximum number of iterations reached without convergence")
            convergence <- 1
            break
        }
        iteration <- iteration + 1

        # Self-adjusting parameter control scheme
        i <- runif(NP) <= tau_F
        # Combine jitter with dither
        Ftrial[, i] <- if (use_jitter)
            (1 + jitter_factor*runif(d, -0.5, 0.5)) %o% runif(sum(i), Fl, Fu)
        else runif(sum(i), Fl, Fu)

        i <- runif(NP) <= tau_CR
        CRtrial[i] <- runif(sum(i), CRl, CRu)

        i <- runif(NP) <= tau_pF
        pFtrial[i] <- runif(sum(i))

        # Equalize the mean lifetime of all vectors
        # Price, KV, Storn, RM, and Lampinen, JA (2005)
        # Differential Evolution: A Practical Approach to
        # Global Optimization. Springer, p 284
        for (i in ((iteration + pop_index) %% NP) + 1) { # Start loop
                                                         # through population
            # DE/rand/1/either-or/bin
            X_i <- pop[, i]
            # Randomly pick 3 vectors all different from target vector
            r <- sample(pop_index[-i], 3)
            X_base <- pop[, r[1L]]
            X_r1 <- pop[, r[2L]]
            X_r2 <- pop[, r[3L]]

            trial[, i] <- handle_bounds(perform_reproduction(), X_base)
        }

        eval(child)

        x_best_ind <- which_best(fpop)

        converge <- eval(conv)
        if (trace && (iteration %% triter == 0))
            cat(iteration, ":", "<", converge, ">", "(", fpop[x_best_ind], ")",
                pop[, x_best_ind],
                if (has_constraints)
                    paste("{", which(hpop[, x_best_ind] > 0), "}"),
                fill = TRUE)
    }

    res <- list(par = pop[, x_best_ind],
                value = fpop[x_best_ind],
                iter = iteration,
                convergence = convergence)
    if (details) {
        res$poppar <- pop
        res$popcost <- fpop
    }
    res
}
