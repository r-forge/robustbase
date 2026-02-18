\name{SPJDEoptim}
\alias{SPJDEoptim}
\title{
  Simple Parallel Version of \code{\link[=JDEoptim]{JDEoptim()}}
}
\description{
  Parallel evaluation of the objective function and/or constraints for each
  candidate solution in the population.
}
\usage{
SPJDEoptim(lower, upper, fn,
           constr = NULL, meq = 0, eps = 1e-5,
           sequential_eval = c("none", "constraints",
                               "objective", "both"),
           ..., further_args = list(),
           NP = 10*length(lower),
           Fl = 0.1, Fu = 1, CRl = 0, CRu = 1.1,
           tau_F = 0.1, tau_CR = 0.1, tau_pF = 0.1,
           jitter_factor = 0.001,
           tol = 1e-15, maxiter = 2000*length(lower),
           fnscale = 1, compare_to = c("median", "max"),
           add_to_init_pop = NULL,
           trace = FALSE, triter = 1,
           details = FALSE)
}
\arguments{
  \item{lower, upper}{numeric vectors, the lower and upper bounds of the
    search space (\emph{box constraints}); must be finite
    (\code{\link{is.finite}}).}

  \item{fn}{a \code{\link{function}} to be \strong{minimized} that takes a
    numeric vector \eqn{X_i} as first argument and returns the value of the
    objective.}

  \item{constr}{a vector \code{\link{function}} specifying the
    \strong{left-hand side} of equality constraints defined to equal zero
    (\eqn{h_j(X_i) = 0,\; j = 1,\ldots,\mathrm{meq}}),
    followed by inequality constraints defined as lesser than or equal to zero
    (\eqn{g_j(X_i) \le 0,\; j = \mathrm{meq}+1,\ldots}).  This function takes
    \eqn{X_i} as its first argument and returns a numeric vector with the same
    length of the total number of constraints.  It defaults to \code{NULL},
    which means that \strong{bound-constrained} minimization is used.}

  \item{meq}{an integer, the first \code{meq} constraints are
    \emph{equality} constraints whereas the remaining ones are
    \emph{inequality} constraints.  Defaults to \code{0}
    (inequality constraints only).}

  \item{eps}{the maximal admissible constraint violation for
    equality constraints.  A numeric vector of small positive tolerance values
    with length \code{meq} used in the transformation of equalities into
    inequalities of the form \eqn{|h_j(X_i)| - \epsilon \le 0}.  A scalar value
    is expanded to apply to all equality constraints.  Default is \code{1e-5}.}

  \item{sequential_eval}{a character string that indicates when either the
    objective function (if \code{"objective"}) or the constraints
    (if \code{"constraints"}) should be evaluated in a serial instead of
    parallel manner whenever their computational cost is low.  If \code{"both"}
    the computation is purely sequential, being intended for debugging and
    testing.  Defaults to \code{"none"}.}

  \item{...}{[\emph{to be used only when run in parallel mode}] named arguments
    specifying \strong{locally-defined} functions and other objects required
    but not defined in \code{fn} and/or \code{constr}.}

  \item{further_args}{a list of any further arguments passed to
    \code{fn} and \code{constr}.  Defaults to an empty list.}

  \item{NP}{an integer, the population size.  Defaults to \code{100}.}

  \item{Fl}{a numeric, the minimum value that the
    \emph{scaling factor} \code{F} could take.  It defaults to \code{0.1}.}

  \item{Fu}{a numeric, the maximum value that the
    \emph{scaling factor} \code{F} could take.  It defaults to \code{1}.}

  \item{CRl}{a numeric, the minimum value to be used for the
    \emph{crossover constant} \code{CR}.  It defaults to \code{0}.}

  \item{CRu}{a numeric, the maximum value to be used for the
    \emph{crossover constant} \code{CR}.  It defaults to \code{1.1}.}

  \item{tau_F}{a numeric, the probability that the
    \emph{scaling factor} \code{F} is updated.  Defaults to \code{0.1}.}

  \item{tau_CR}{a numeric, the probability that the
    \emph{crossover constant} \code{CR} is updated.  Defaults to \code{0.1}.}

  \item{tau_pF}{a numeric, the probability that the
    \emph{mutation probability} \eqn{p_F}{pF} in the mutation strategy
    DE/rand/1/either-or is updated.  Defaults to \code{0.1}.}

  \item{jitter_factor}{a numeric, the tuning constant for \emph{jitter}.
    If \code{NULL} only \emph{dither} is used.  Default is \code{0.001}.}

  \item{tol}{a numeric, the tolerance for the stopping criterion.  Default is
    \code{1e-15}.}

  \item{maxiter}{an integer, the maximum number of iterations allowed.
    Default is \code{2000*length(lower)}.}

  \item{fnscale}{a numeric, the typical magnitude of \code{fn} used in the
    \emph{stopping criterion}.  Defaults to \code{1}.  See \sQuote{Details} in
    \code{\link[=JDEoptim]{JDEoptim()}}.}

  \item{compare_to}{a character string controlling which function should be
    applied when evaluating the \emph{stopping criterion}.  It defaults to
    \code{"median"}.  See \sQuote{Details} in
    \code{\link[=JDEoptim]{JDEoptim()}}.}

  \item{add_to_init_pop}{numeric vector of length \code{length(lower)} or
    column-wise \code{\link{matrix}} with \code{length(lower)} rows specifying
    initial candidate solutions which are appended to the randomly generated
    initial population.  Default is \code{NULL}.}

  \item{trace}{logical indicating whether or not to monitor the
    iteration process.  Default is \code{FALSE}.}

  \item{triter}{an integer, trace output is printed at every
    \code{triter} iterations.  Default is \code{1}.}

  \item{details}{logical indicating if the output will contain the solutions
    in the final population and their respective \code{fn} values.
    Defaults to \code{FALSE}.}
}
\details{
  The optimizer is intended for dealing with computationally expensive
  objective functions and/or constraints (see, \emph{e.g.},
  Storn and Price, 2025).  The algorithm has been parallelized simply by
  distributing the evaluation of the candidate solutions in the population over
  the available processors/cores with
  \code{\link[=mirai_map]{mirai_map()}} from package \pkg{mirai}.

  It is important to note this approach and its implementation
  \strong{is not equivalent} to the purely sequential version
  \code{JDEoptim()}, the major reason being the fact that an asynchronous
  update of the current population is not possible anymore.  This in turn almost
  always has a \strong{deleterious effect} on the efficiency of DE.  For this
  reason although it can run \strong{without \pkg{mirai} installed}, for example
  using \code{sequential_eval = "both"}, the use of
  \code{\link[=JDEoptim]{JDEoptim()}} is preferred over \code{SPJDEoptim()}
  in this case.

  Note that all support for the parallel computation is provided by package
  \CRANpkg{mirai} and controlled directly by end-users; in particular,
  resources are declared with \code{\link[=daemons]{mirai::daemons()}}
  \strong{before invoking} \code{SPJDEoptim()}.
}
\value{
  A list with the following components:
  \item{par}{the best solution \eqn{X_i} found.}

  \item{value}{the value of \eqn{\mathrm{fn}(X_i)} corresponding to \code{par}.}

  \item{iter}{the number of iterations used.}

  \item{convergence}{an integer code; \code{0} indicates successful completion,
    whereas \code{1} indicates that the iteration limit \code{maxiter} has been
    reached.}

  and if \code{details = TRUE}:
  \item{poppar}{a \code{\link{matrix}} whose columns are the solutions remaining
    in the final population.}

  \item{popcost}{the values of \eqn{\mathrm{fn}(X_i)} for the
    corresponding columns of \code{poppar}.}
}
\references{
  \url{https://mirai.r-lib.org}

  Storn, Rainer and Price, Kenneth V. (2025)
  Multi-child DE --- a massively parallel differential evolution algorithm;
  in \emph{2025 IEEE Symposium on Computational Intelligence on
  Engineering/Cyber Physical Systems (CIES)}, Trondheim, Norway.
  \doi{10.1109/CIES64955.2025.11007629}.
}
\author{
  Eduardo L. T. Conceicao \email{mail@eduardoconceicao.org}
}
\note{
  Additional arguments are passed to \code{fn} and \code{constr} through
  \code{further_args} instead of \code{...} such as in
  \code{\link[=JDEoptim]{JDEoptim()}}.
}
\examples{
\donttest{
if (requireNamespace("mirai", quietly = TRUE)) {
# NOTE: Examples were excluded from testing
#       to reduce package check time.

# stays within 2-core limit
mirai::daemons(1, dispatcher = FALSE)
set.seed(42)

# Warning: the examples are run using a very small number of
# iterations to decrease execution time. Because they are very
# cheap to compute, parallel evaluation is actually detrimental.

#   Shekel 5 problem
#
#   0 <= xj <= 10, j = {1, 2, 3, 4}
#   There are five local minima and the global minimizer is
#   located at x* = (4.00, 4.00, 4.00, 4.00) with
#   f(x*) ~~ -10.1499.
#
#   Source:
#     Ali, M. Montaz, Khompatraporn, Charoenchai, and
#     Zabinsky, Zelda B. (2005).
#     A numerical evaluation of several stochastic algorithms
#     on selected continuous global optimization test problems.
#     Journal of Global Optimization 31, 635-672.
#     https://doi.org/10.1007/s10898-004-9972-2
a_ij <- matrix(c(4, 4, 4, 4,
                 1, 1, 1, 1,
                 8, 8, 8, 8,
                 6, 6, 6, 6,
                 3, 7, 3, 7), ncol = 4, byrow = TRUE)
c_i <- c(0.1, 0.2, 0.2, 0.4, 0.4)
a_ji <- t(a_ij)
shekel5 <- function(x) {
    -sum(1/( colSums((x - a_ji)^2) + c_i ))
}

SPJDEoptim(rep(0, 4), rep(10, 4), shekel5,
           a_ji = a_ji, c_i = c_i,
           maxiter = 10, trace = TRUE)

#   Sinusoidal problem
#
#   0 <= xi <= 180, i = {1, 2, ..., n}
#   The variable x is in degrees. Parameter A affects the
#   amplitude of the global optimum; B affects the peridiocity
#   and hence the number of local minima; z shifts the location
#   of the global minimum. The location of the global solution
#   is at x* = (90+z, 90+z, ..., 90+z) with the global optimum
#   value of f(x*) = -(A+1).
#
#   Source:
#     Ali, M. Montaz, Khompatraporn, Charoenchai, and
#     Zabinsky, Zelda B. (2005).
#     A numerical evaluation of several stochastic algorithms
#     on selected continuous global optimization test problems.
#     Journal of Global Optimization 31, 635-672.
#     https://doi.org/10.1007/s10898-004-9972-2
sinusoidal <- function(x, A, B, z) {
    var_in_rad <- (x - z)*pi/180
    -( A*prod(sin(var_in_rad)) + prod(sin(B*var_in_rad)) )
}

SPJDEoptim(rep(0, 10), rep(180, 10), sinusoidal,
           further_args = list(A = 2.5, B = 5, z = 30),
           maxiter = 10, trace = TRUE)

mirai::daemons(0) # reset
Sys.sleep(1)
}
}
}
\concept{parallel computation}
