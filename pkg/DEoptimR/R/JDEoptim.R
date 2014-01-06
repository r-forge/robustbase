##' Nonlinear Constrained Optimization via Differential Evolution
##'
##' An implementation of the jDE variant of the Differential Evolution
##' stochastic algorithm for global optimization of nonlinear programming
##' problems.
##'
##' The setting of the \emph{control parameters} of standard Differential
##' Evolution (DE) is crucial for the algorithm's performance. Unfortunately,
##' when the generally recommended values for these parameters (see,
##' \emph{e.g.}, Storn and Price, 1997) are unsuitable for use, their
##' determination is often difficult and time consuming. The jDE algorithm
##' proposed in Brest \emph{et al.} (2006) employs a simple self-adaptive
##' scheme to perform the automatic setting of control parameters scale factor
##' \code{F} and crossover rate \code{CR}.
##'
##' This implementation differs from the original description, most notably in
##' the use of the \emph{DE/rand/1/either-or} mutation strategy (Price \emph{et
##' al.}, 2005), combination of \emph{jitter with dither} (Storn 2008), and its
##' use of only a \emph{single population} (Babu and Angira 2006) instead of
##' separate current and child populations as in classical DE.
##'
##' Constraint handling is done using the approach described in Zhang and
##' Rangaiah (2012).
##'
##' DE is easily extended to deal with \emph{mixed integer nonlinear
##' programming} problems using a small variation of the technique presented by
##' Lampinen and Zelinka (1999). Integer values are obtained by means of the
##' \code{floor()} function for the evaluation of the objective function. This
##' is because DE itself works with continuous variables. Additionally, each
##' upper bound of the integer variables should be added by \code{1}.
##'
##' The algorithm is stopped if
##' \deqn{\frac{\mathrm{FUN}\{[\mathrm{fn}(x_1),\ldots,\mathrm{fn}(x_\mathrm{npop})]\}
##' - \mathrm{fn}(x_\mathrm{best})}{\mathrm{fnscale}} \le \mathrm{tol}}{% (
##' FUN{ [fn(x[1]),\ldots,fn(x[npop])] } - fn(x[best]) )/fnscale <= tol,} where
##' the \dQuote{best} individual \eqn{x_\mathrm{best}}{x[best]} is the
##' \emph{feasible} solution with the lowest objective function value in the
##' population and the total number of elements in the population, \code{npop},
##' is \code{NP+NCOL(add_to_init_pop)}.  This is a variant of the \emph{Diff}
##' criterion studied by Zielinski and Laur (2008), which was found to yield
##' the best results.
##'
##' @param lower,upper numeric vectors of \emph{lower} or \emph{upper} bounds,
##' respectively, for the parameters to be optimized over.  Must be finite
##' (\code{\link{is.finite}}) as they bound the hyper rectangle of the initial
##' random population.
##' @param fn (nonlinear) objective \code{\link{function}} to be
##' \emph{minimized}.  It takes as first argument the vector of parameters over
##' which minimization is to take place.  It must return the value of the
##' function at that point.
##' @param constr an optional \code{\link{function}} for specifying the
##' nonlinear constraints under which we want to minimize \code{fn}.  They
##' should be given in the form \eqn{h_i(x) = 0, g_i(x) \le 0}.  This function
##' takes the vector of parameters as its first argument and returns a real
##' vector with the length of the total number of constraints.  It defaults to
##' \code{NULL}, meaning that \emph{bound-constrained} minimization is used.
##' @param meq an optional positive integer specifying that the first
##' \code{meq} constraints are treated as \emph{equality} constraints, all the
##' remaining as \emph{inequality} constraints. Defaults to \code{0}
##' (inequality constraints only).
##' @param eps an optional real vector of small positive tolerance values with
##' length \code{meq} used in the transformation of equalities into
##' inequalities of the form \eqn{|h_i(x)| - \epsilon \le 0}. A scalar value is
##' expanded to apply to all equality constraints. Default is \code{1e-5}.
##' @param NP an optional positive integer giving the number of candidate
##' solutions in the randomly distributed initial population. Defaults to
##' \code{10*length(lower)}.
##' @param Fl an optional scalar which represents the minimum value that the
##' \emph{scaling factor} \code{F} could take. Default is \code{0.1}, which is
##' almost always satisfactory.
##' @param Fu an optional scalar which represents the maximum value that the
##' \emph{scaling factor} \code{F} could take. Default is \code{1}, which is
##' almost always satisfactory.
##' @param tau1 an optional scalar which represents a probability in the
##' mutation strategy DE/rand/1/either-or. Defaults to \code{0.1}.
##' @param tau2 an optional scalar which represents the probability that the
##' \emph{scaling factor} \code{F} is updated. Defaults to \code{0.1}, which is
##' almost always satisfactory.
##' @param tau3 an optional constant value which represents the probability
##' that the \emph{crossover probability} \code{CR} is updated. Defaults to
##' \code{0.1}, which is almost always satisfactory.
##' @param jitter_factor an optional tuning constant for \emph{jitter}. If
##' \code{NULL} only \emph{dither} is used. Defaults to \code{0.001}.
##' @param tol an optional positive scalar giving the tolerance for the
##' stopping criterion. Default is \code{1e-15}.
##' @param maxiter an optional positive integer specifying the maximum number
##' of iterations that may be performed before the algorithm is halted.
##' Defaults to \code{200*length(lower)}.
##' @param fnscale an optional positive scalar specifying the typical magnitude
##' of \code{fn}. It is used only in the \emph{stopping criterion}. Defaults to
##' \code{1}. See \sQuote{Details}.
##' @param FUN an optional character string controlling which function should
##' be applied to the \code{fn} values of the candidate solutions in a
##' generation to be compared with the so-far best one when evaluating the
##' \emph{stopping criterion}. If \dQuote{\code{median}} the \code{median}
##' function is used; else, if \dQuote{\code{max}} the \code{max} function is
##' used. It defaults to \dQuote{\code{median}}. See \sQuote{Details}.
##' @param add_to_init_pop an optional real vector of length
##' \code{length(lower)} or matrix with \code{length(lower)} rows specifying
##' initial values of the parameters to be optimized which are appended to the
##' randomly generated initial population. It defaults to \code{NULL}.
##' @param trace an optional logical value indicating if a trace of the
##' iteration progress should be printed. Default is \code{FALSE}.
##' @param triter an optional positive integer that controls the frequency of
##' tracing when \code{trace = TRUE}. Default is \code{triter = 1}, which means
##' that \code{iteration : < value of stopping test > ( value of best solution
##' ) best solution { index of violated constraints }} is printed at every
##' iteration.
##' @param details an optional logical value. If \code{TRUE} the output will
##' contain the parameters in the final population and their respective
##' \code{fn} values. Defaults to \code{FALSE}.
##' @param ... optional additional arguments passed to \code{fn()} \emph{and}
##' \code{constr()} if that is not \code{NULL}.
##' @return A list with the following components: \item{par}{The best set of
##' parameters found.}
##'
##' \item{value}{The value of \code{fn} corresponding to \code{par}.}
##'
##' \item{iter}{Number of iterations taken by the algorithm.}
##'
##' \item{convergence}{An integer code. \code{0} indicates successful
##' completion. \code{1} indicates that the iteration limit \code{maxiter} has
##' been reached.} and if \code{details = TRUE}: \item{poppar}{Matrix of
##' dimension \code{(length(lower), npop)}, with columns corresponding to the
##' parameter vectors remaining in the population.}
##'
##' \item{popcost}{The values of \code{fn} associated with \code{poppar},
##' vector of length \code{npop}.}
##' @note It is possible to perform a warm start, \emph{i.e.}, starting from
##' the previous run and resume optimization, using \code{NP = 0} and the
##' component \code{poppar} for the \code{add_to_init_pop} argument.
##' @author Eduardo L. T. Conceicao \email{econceicao@@kanguru.pt}
##' @seealso Function \code{\link[DEoptim]{DEoptim}()} in the \pkg{DEoptim}
##' package has many more options than \code{JDEoptim()}, but does not allow
##' constraints in the same flexible manner.
##' @references Babu, B. V. and Angira, R. (2006) Modified differential
##' evolution (MDE) for optimization of non-linear chemical processes.
##' \emph{Computers and Chemical Engineering} \bold{30}, 989--1002.
##'
##' Brest, J., Greiner, S., Boskovic, B., Mernik, M. and Zumer, V. (2006)
##' Self-adapting control parameters in differential evolution: a comparative
##' study on numerical benchmark problems.  \emph{IEEE Transactions on
##' Evolutionary Computation} \bold{10}, 646--657.
##'
##' Lampinen, J. and Zelinka, I. (1999).  Mechanical engineering design
##' optimization by differential evolution; in Corne, D., Dorigo, M. and
##' Glover, F., Eds., \emph{New Ideas in Optimization}.  McGraw-Hill, pp.
##' 127--146.
##'
##' Price, K. V., Storn, R. M. and Lampinen, J. A. (2005) \emph{Differential
##' Evolution: A practical approach to global optimization}.  Springer, Berlin,
##' pp. 117--118.
##'
##' Storn, R. (2008) Differential evolution research --- trends and open
##' questions; in Chakraborty, U. K., Ed., \emph{Advances in differential
##' evolution}.  SCI 143, Springer-Verlag, Berlin, pp. 11--12.
##'
##' Storn, R. and Price, K. (1997) Differential evolution - a simple and
##' efficient heuristic for global optimization over continuous spaces.
##' \emph{Journal of Global Optimization} \bold{11}, 341--359.
##'
##' Zhang, H. and Rangaiah, G. P. (2012) An efficient constraint handling
##' method with integrated differential evolution for numerical and engineering
##' optimization.  \emph{Computers and Chemical Engineering} \bold{37}, 74--88.
##'
##' Zielinski, K. and Laur, R. (2008) Stopping criteria for differential
##' evolution in constrained single-objective optimization; in Chakraborty, U.
##' K., Ed., \emph{Advances in differential evolution}.  SCI 143,
##' Springer-Verlag, Berlin, pp. 111--138.
##' @examples
##'
##' # Use a preset seed so test values are reproducible.
##' set.seed(1234)
##'
##' # Bound-constrained optimization
##'
##' #   Griewank function
##' #
##' #   -600 <= xi <= 600, i = {1, 2, ..., n}
##' #   The function has a global minimum located at
##' #   x* = (0, 0, ..., 0) with f(x*) = 0. Number of local minima
##' #   for arbitrary n is unknown, but in the two dimensional case
##' #   there are some 500 local minima.
##' #
##' #   Source:
##' #     Ali, M. Montaz, Khompatraporn, Charoenchai, and
##' #     Zabinsky, Zelda B. (2005).
##' #     A numerical evaluation of several stochastic algorithms
##' #     on selected continuous global optimization test problems.
##' #     Journal of Global Optimization 31, 635-672.
##' griewank <- function(x) {
##'     1 + crossprod(x)/4000 - prod( cos(x/sqrt(seq_along(x))) )
##' }
##'
##' JDEoptim(rep(-600, 10), rep(600, 10), griewank,
##'          tol = 1e-7, trace = TRUE, triter = 50)
##'
##' # Nonlinear constrained optimization
##'
##' #   0 <= x1 <= 34, 0 <= x2 <= 17, 100 <= x3 <= 300
##' #   The global optimum is
##' #   (x1, x2, x3; f) = (0, 16.666667, 100; 189.311627).
##' #
##' #   Source:
##' #     Westerberg, Arthur W., and Shah, Jigar V. (1978).
##' #     Assuring a global optimum by the use of an upper bound
##' #     on the lower (dual) bound.
##' #     Computers and Chemical Engineering 2, 83-92.
##' fcn <-
##'     list(obj = function(x) {
##'              35*x[1]^0.6 + 35*x[2]^0.6
##'          },
##'          eq = 2,
##'          con = function(x) {
##'              x1 <- x[1]; x3 <- x[3]
##'              c(600*x1 - 50*x3 - x1*x3 + 5000,
##'                600*x[2] + 50*x3 - 15000)
##'          })
##'
##' JDEoptim(c(0, 0, 100), c(34, 17, 300),
##'          fn = fcn$obj, constr = fcn$con, meq = fcn$eq,
##'          tol = 1e-7, trace = TRUE, triter = 50)
##'
##' #   Designing a pressure vessel
##' #   Case A: all variables are treated as continuous
##' #
##' #   1.1 <= x1 <= 12.5*, 0.6 <= x2 <= 12.5*,
##' #   0.0 <= x3 <= 240.0*, 0.0 <= x4 <= 240.0
##' #   Roughly guessed*
##' #   The global optimum is (x1, x2, x3, x4; f) =
##' #   (1.100000, 0.600000, 56.99482, 51.00125; 7019.031).
##' #
##' #   Source:
##' #     Lampinen, Jouni, and Zelinka, Ivan (1999).
##' #     Mechanical engineering design optimization
##' #     by differential evolution.
##' #     In: David Corne, Marco Dorigo and Fred Glover (Editors),
##' #     New Ideas in Optimization, McGraw-Hill, pp 127-146
##' pressure_vessel_A <-
##'     list(obj = function(x) {
##'              x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]
##'              0.6224*x1*x3*x4 + 1.7781*x2*x3^2 +
##'              3.1611*x1^2*x4 + 19.84*x1^2*x3
##'          },
##'          con = function(x) {
##'              x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]
##'              c(0.0193*x3 - x1,
##'                0.00954*x3 - x2,
##'                750.0*1728.0 - pi*x3^2*x4 - 4/3*pi*x3^3)
##'          })
##'
##' JDEoptim(c( 1.1,  0.6,   0.0,   0.0),
##'          c(12.5, 12.5, 240.0, 240.0),
##'          fn = pressure_vessel_A$obj,
##'          constr = pressure_vessel_A$con,
##'          tol = 1e-7, trace = TRUE, triter = 50)
##'
##' # Mixed integer nonlinear programming
##'
##' #   Designing a pressure vessel
##' #   Case B: solved according to the original problem statements
##' #           steel plate available in thicknesses multiple
##' #           of 0.0625 inch
##' #
##' #   wall thickness of the
##' #   shell 1.1 [18*0.0625] <= x1 <= 12.5 [200*0.0625]
##' #   heads 0.6 [10*0.0625] <= x2 <= 12.5 [200*0.0625]
##' #         0.0 <= x3 <= 240.0, 0.0 <= x4 <= 240.0
##' #   The global optimum is (x1, x2, x3, x4; f) =
##' #   (1.125 [18*0.0625], 0.625 [10*0.0625],
##' #    58.29016, 43.69266; 7197.729).
##' pressure_vessel_B <-
##'     list(obj = function(x) {
##'              x1 <- floor(x[1])*0.0625
##'              x2 <- floor(x[2])*0.0625
##'              x3 <- x[3]; x4 <- x[4]
##'              0.6224*x1*x3*x4 + 1.7781*x2*x3^2 +
##'              3.1611*x1^2*x4 + 19.84*x1^2*x3
##'          },
##'          con = function(x) {
##'              x1 <- floor(x[1])*0.0625
##'              x2 <- floor(x[2])*0.0625
##'              x3 <- x[3]; x4 <- x[4]
##'              c(0.0193*x3 - x1,
##'                0.00954*x3 - x2,
##'                750.0*1728.0 - pi*x3^2*x4 - 4/3*pi*x3^3)
##'          })
##'
##' JDEoptim(c( 18,    10,     0.0,   0.0),
##'          c(200+1, 200+1, 240.0, 240.0),
##'          fn = pressure_vessel_B$obj,
##'          constr = pressure_vessel_B$con,
##'          tol = 1e-7, trace = TRUE, triter = 50)
##'
JDEoptim <-
    function(lower, upper, fn, constr = NULL, meq = 0, eps = 1e-5,
             NP = 10*d, Fl = 0.1, Fu = 1, tau1 = 0.1, tau2 = 0.1, tau3 = 0.1,
             jitter_factor = 0.001,
             tol = 1e-15, maxiter = 200*d, fnscale = 1,
             FUN = c("median", "max"),
             add_to_init_pop = NULL, trace = FALSE, triter = 1,
             details = FALSE, ...)

#   Copyright 2013, Eduardo L. T. Conceicao
#   Available under the GPL (>= 2)

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
        trial <- if (runif(1) <= pFtrial)
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
                if (TAVtrial <= TAVpop[i]) { # trial and target are both
                    pop[, i] <- trial        # unfeasible, the one with smaller
                    hpop[, i] <- htrial      # constraint violation is chosen
                    F[, i] <- Ftrial         # or trial vector when both are
                    CR[i] <- CRtrial         # solutions of equal quality
                    pF[i] <- pFtrial
                    TAVpop[i] <- TAVtrial
                }
            } else if (TAVpop[i] > mu) { # trial is feasible and target is not
                pop[, i] <- trial
                fpop[i] <- fn1(trial)
                hpop[, i] <- htrial
                F[, i] <- Ftrial
                CR[i] <- CRtrial
                pF[i] <- pFtrial
                TAVpop[i] <- TAVtrial
                FF <- sum(TAVpop <= mu)/NP
                mu <- mu*(1 - FF/NP)
            } else {                     # between two feasible solutions, the
                ftrial <- fn1(trial)     # one with better objective function
                if (ftrial <= fpop[i]) { # value is chosen
                    pop[, i] <- trial    # or trial vector when both are
                    fpop[i] <- ftrial    # solutions of equal quality
                    hpop[, i] <- htrial
                    F[, i] <- Ftrial
                    CR[i] <- CRtrial
                    pF[i] <- pFtrial
                    TAVpop[i] <- TAVtrial
                    FF <- sum(TAVpop <= mu)/NP
                    mu <- mu*(1 - FF/NP)
                }
            }
        })
    }

    which.best <-
        if (!is.null(constr))
            function(x) {
                ind <- TAVpop <= mu
                if (all(ind))
                    which.min(x)
                else if (any(ind))
                    which(ind)[which.min(x[ind])]
                else which.min(TAVpop)
            }
        else which.min

    # Check input parameters
    FUN <- match.arg(FUN)
    d <- length(lower)
    if (length(upper) != d)
        stop("'lower' must have same length as 'upper'")
    stopifnot(is.numeric(lower), is.numeric(upper),
	      is.finite(lower), is.finite(upper), lower <= upper,
	      length(fnscale) == 1, is.finite(fnscale), fnscale > 0,
	      is.function(fn))
    if (!is.null(constr)) {
        stopifnot(is.function(constr))
        stopifnot(length(meq) == 1, meq == as.integer(meq), meq >= 0,
                  is.numeric(eps), is.finite(eps), eps > 0)
        if (length(eps) == 1)
            eps <- rep_len(eps, meq)
        else if (length(eps) != meq)
            stop("eps must be either of length meq, or length 1")
    }
    stopifnot(length(NP) == 1, NP == as.integer(NP),
	      length(Fl) == 1, is.numeric(Fl),
	      length(Fu) == 1, is.numeric(Fu),
	      Fl <= Fu)
    stopifnot(length(tau1) == 1, is.numeric(tau1), 0 <= tau1, tau1 <= 1,
	      length(tau2) == 1, is.numeric(tau2), 0 <= tau2, tau2 <= 1,
	      length(tau3) == 1, is.numeric(tau3), 0 <= tau3, tau3 <= 1)
    if (!is.null(jitter_factor))
        stopifnot(length(jitter_factor) == 1, is.numeric(jitter_factor))
    stopifnot(length(tol) == 1, is.numeric(tol),
              length(maxiter) == 1, maxiter == as.integer(maxiter),
              length(triter) == 1, triter == as.integer(triter))
    if (!is.null(add_to_init_pop))
        stopifnot(NROW(add_to_init_pop) == d,
                  is.numeric(add_to_init_pop),
                  add_to_init_pop >= lower,
                  add_to_init_pop <= upper)

    # Initialization:
    fn1 <- function(par) fn(par, ...)

    if (!is.null(constr))
        constr1 <-
            if (meq > 0) {
                eqI <- 1:meq
                function(par) {
                    h <- constr(par, ...)
                    h[eqI] <- abs(h[eqI]) - eps
                    h
                }
            } else function(par) constr(par, ...)

    use.jitter <- !is.null(jitter_factor)

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
    stopifnot(NP >= 4)
    F <- if (use.jitter)
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
        if (iteration >= maxiter) {
            warning("maximum number of iterations reached without convergence")
            convergence <- 1
            break
        }
        iteration <- iteration + 1

        for (i in popIndex) { # Start loop through population
            # Equalize the mean lifetime of all vectors
            # Price, KV, Storn, RM, and Lampinen, JA (2005)
            # Differential Evolution. Springer, p 284
            i <- ((iteration + i) %% NP) + 1

            # Fi update
            # Combine jitter with dither
            # Storn, Rainer (2008).
            # Differential evolution research - trends and open questions.
            # In: U. K. Chakraborty (Ed.), Advances in Differential Evolution,
            # SCI 143, Springer-Verlag, pp 11-12
            Ftrial <- if (runif(1) <= tau1) {
                if (use.jitter)
                    runif(1, Fl, Fu) * (1 + jitter_factor*runif(d, -0.5, 0.5))
                else runif(1, Fl, Fu)
            } else F[, i]
            # CRi update
            CRtrial <- if (runif(1) <= tau2) runif(1) else CR[i]
            # pFi update
            pFtrial <- if (runif(1) <= tau3)
                runif(1)
            else pF[i]

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
                iter = iteration,
                convergence = convergence)
    if (details) {
        res$poppar <- pop
        res$popcost <- fpop
    }
    res
}

## Not exported, and only used because CRAN checks must be faster
doExtras <- function() {
    interactive() || nzchar(Sys.getenv("R_JDEoptim_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}
