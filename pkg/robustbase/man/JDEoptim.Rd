\name{JDEoptim}
\alias{JDEoptim}
\title{Nonlinear Constrained Optimization via Differential Evolution}
\description{
  An implementation of the jDE variant of the Differential Evolution stochastic algorithm for global optimization of nonlinear programming problems.
}
\usage{
  JDEoptim(lower, upper, fn,
    constr = NULL, meq = 0, eps = 1e-05,
    NP = 10*d, Fl = 0.1, Fu = 1,
    tau1 = 0.1, tau2 = 0.1, tau3 = 0.1,
    jitter_factor = 0.001,
    tol = 1e-15, maxiter = 200*d, fnscale = 1,
    FUN = c("median", "max"),
    add_to_init_pop = NULL, trace = FALSE, triter = 1,
    details = FALSE, ...)
}
\arguments{
  \item{lower}{A real vector of \emph{lower} bounds for the parameters to be optimized over.}
  \item{upper}{A real vector of \emph{upper} bounds for the parameters to be optimized over.}

  \item{fn}{A nonlinear objective \code{\link{function}} to be \emph{minimized} that takes as first argument the vector of parameters over which minimization is to take place. It should return a scalar result that is the value of the function at that point.}

  \item{constr}{An optional \code{\link{function}} for specifying the nonlinear constraints under which we want to minimize \code{fn}. They should be given in the form \eqn{h_i(x) = 0, g_i(x) <= 0}. This function takes the vector of parameters as its first argument and returns a real vector with the length of the total number of constraints. It defaults to \code{NULL}, meaning that \emph{bound-constrained} minimization is used.}

  \item{meq}{An optional positive integer specifying that the first \code{meq} constraints are treated as \emph{equality} constraints, all the remaining as \emph{inequality} constraints. Defaults to \code{0} (inequality constraints only).}

  \item{eps}{An optional real vector of small positive tolerance values with length \code{meq} which transforms equalities into \eqn{|h_i(x)| - \epsilon <= 0} inequalities. A scalar value is expanded to apply to all equalities. Default is \code{1e-5}.}

  \item{NP}{An optional positive integer giving the number of candidate solutions in the randomly distributed initial population. Defaults to \code{10*length(lower)}.}

  \item{Fl}{An optional scalar which represents the minimum value that the \emph{scaling factor} \code{F} could take. The default is almost always satisfactory. Default is \code{0.1}.}

  \item{Fu}{An optional scalar which represents the maximum value that the \emph{scaling factor} \code{F} could take. The default is almost always satisfactory. Default is \code{1}.}

  \item{tau1}{An optional scalar which represents a probability in the mutation strategy DE/rand/1/either-or. Defaults to \code{0.1}.}

  \item{tau2}{An optional scalar which represents the probability that the \emph{scaling factor} \code{F} is updated. The default is almost always satisfactory. Defaults to \code{0.1}.}

  \item{tau3}{An optional constant value which represents the probability that the \emph{crossover probability} \code{CR} is updated. The default is almost always satisfactory. Defaults to \code{0.1}.}

  \item{jitter_factor}{An optional tuning constant for \emph{jitter}. If \code{NULL} only \emph{dither} is used. Defaults to \code{0.001}.}

  \item{tol}{An optional positive scalar giving the tolerance for the stopping criterion. Default is \code{1e-15}.}

  \item{maxiter}{An optional positive integer specifying the maximum number of iterations that may be performed before the algorithm is halted. Defaults to \code{200*length(lower)}.}

  \item{fnscale}{An optional positive scalar specifying the typical magnitude of \code{fn}. It is used only in the stopping criterion given below in \sQuote{Details}. Defaults to \code{1}.}

  \item{FUN}{An optional character string controlling which function should be applied to the population of candidate solutions \code{fn} values to be compared with the so-far best one \code{min(fn(.))}. If \dQuote{\code{median}} the \code{median} function is used; else, if \dQuote{\code{max}} the \code{max} function is used. It defaults to \dQuote{\code{median}}.}

  \item{add_to_init_pop}{An optional real vector of length \code{length(lower)} or matrix with \code{length(lower)} rows specifying initial values of the parameters to be optimized which are appended to the randomly generated initial population. It defaults to \code{NULL}.}

  \item{trace}{An optional logical value indicating if a trace of the iteration progress should be printed. Default is \code{FALSE}.}

  \item{triter}{An optional positive integer that controls the frequency of tracing when \code{trace = TRUE}. Default is \code{triter = 1}, which means that \code{iteration : < value of stopping test > ( value of best solution ) best solution { index of violated constraints }} is printed at every iteration.}

  \item{details}{An optional logical value. If \code{TRUE} the output will contain the parameters in the final population and their respective \code{fn} values. Defaults to \code{FALSE}.}

  \item{...}{Optional additional arguments passed to \code{fn()}
  \emph{and} \code{constr()} if that is not \code{NULL}.}
}
\details{
  The setting of the \emph{control parameters} of standard Differential Evolution (DE) is crucial for the algorithm's performance. Unfortunately, when the generally recommended values for these parameters (see, \emph{e.g.}, Storn and Price, 1997) are unsuitable for use, their determination is often difficult and time consuming. The jDE algorithm proposed in Brest \emph{et al.} (2006) employs a simple self-adaptive scheme to perform the automatic setting of control parameters scale factor \code{F} and crossover rate \code{CR}.

  This implementation differs from the original description, most notably in the use of the \emph{DE/rand/1/either-or} mutation strategy (Price \emph{et al.}, 2005), combination of \emph{jitter with dither} (Storn 2008), and its use of only a \emph{single population} (Babu and Angira 2006) instead of two as in classical DE.

  Constraint handling is done using the approach described in Zhang and Rangaiah (2012).

  MORE HERE
}
\value{
  A list with the following components:
  \item{par}{The best set of parameters found.}

  \item{value}{The value of \code{fn} corresponding to \code{par}.}

  \item{iter}{Number of iterations taken by the algorithm.}

  \item{convergence}{An integer code. \code{0} indicates successful completion. \code{1} indicates that the iteration limit \code{maxiter} has been reached.}
  and if \code{details = TRUE}:
  \item{poppar}{Matrix of dimension \code{(length(lower), npop)}, with columns corresponding to the parameter vectors remaining in the population. The total number of elements in the population, \code{npop}, is \code{NP + NCOL(add_to_init_pop)}.}

  \item{popcost}{The values of \code{fn} associated with \code{poppar}, vector of length \code{npop}.}
}
\note{
  It is possible to perform a warm start, \emph{i.e.}, starting from the previous run and resume optimization, using \code{NP = 0} and the component \code{poppar} for the \code{add_to_init_pop} argument.
}
\author{
  Eduardo L. T. Conceicao \email{econceicao@kanguru.pt}
}
\references{
  Babu, B. V. and Angira, R. (2006)
  Modified differential evolution (MDE) for optimization of
  non-linear chemical processes.
  \emph{Computers and Chemical Engineering}, \bold{30}, 989--1002.

  Brest, J., Greiner, S., Boskovic, B., Mernik, M. and Zumer, V. (2006)
  Self-adapting control parameters in differential evolution: a
  comparative study on numerical benchmark problems.
  \emph{IEEE Transactions on Evolutionary Computation}, \bold{10}, 646--657.

  Price, K. V., Storn, R. M. and Lampinen, J. A. (2005)
  \emph{Differential Evolution: A practical approach to global optimization}.
  Springer, Berlin, pp. 117--118.

  Storn, R. (2008)
  Differential evolution research --- trends and open questions;
  in Chakraborty, U. K., Ed.,
  \emph{Advances in differential evolution}.
  SCI 143, Springer-Verlag, Berlin, pp. 11--12.

  Storn, R. and Price, K. (1997)
  Differential evolution - a simple and efficient heuristic for global
  optimization over continuous spaces.
  \emph{Journal of Global Optimization}, \bold{11}, 341--359.

  Zhang, H. and Rangaiah, G. P. (2012)
  An efficient constraint handling method with integrated differential
  evolution for numerical and engineering optimization.
  \emph{Computers and Chemical Engineering}, \bold{37}, 74--88.
}
\seealso{
  Function \code{\link[DEoptim]{DEoptim()}} in the
  \pkg{DEoptim} package has many more options than
  \code{JDEoptim()}, but does not allow constraints in the same
  flexible manner.
}
\examples{
# Use a preset seed so test values are reproducible.
set.seed(1234)

# Bound-constrained optimization

#   Griewank function
#
#   -600 <= xi <= 600, i = {1, 2, ..., n}
#   The global minimum is located at x* = (0, 0, ..., 0)
#   with f(x*) = 0.
griewank <- function(x) {
    1 + crossprod(x)/4000 - prod( cos(x/sqrt(seq_along(x))) )
}

JDEoptim(rep(-600, 10), rep(600, 10), griewank,
         tol = 1e-7, trace = TRUE, triter = 50)

# Nonlinear constrained optimization

#   0 <= x1 <= 34, 0 <= x2 <= 17, 100 <= x3 <= 300
#   The global optimum is
#   (x1, x2, x3; f) = (0, 16.666667, 100; 189.311627).
#
#   Source:
#     Westerberg, Arthur W., and Shah, Jigar V. (1978).
#     Assuring a global optimum by the use of an upper bound
#     on the lower (dual) bound.
#     Computers and Chemical Engineering 2, 83-92.
fcn <-
    list(obj = function(x) {
             35*x[1]^0.6 + 35*x[2]^0.6
         },
         eq = 2,
         con = function(x) {
             x1 <- x[1]; x3 <- x[3]
             c(600*x1 - 50*x3 - x1*x3 + 5000,
               600*x[2] + 50*x3 - 15000)
         })

JDEoptim(c(0, 0, 100), c(34, 17, 300),
         fn = fcn$obj, constr = fcn$con, meq = fcn$eq,
         tol = 1e-7, trace = TRUE, triter = 50)

#   Designing a pressure vessel
#   Case A: all variables are treated as continuous
#
#   1.1 <= x1 <= 12.5*, 0.6 <= x2 <= 12.5*,
#   0.0 <= x3 <= 240.0*, 0.0 <= x4 <= 240.0
#   Roughly guessed*
#   The global optimum is (x1, x2, x3, x4; f) =
#   (1.100000, 0.600000, 56.99482, 51.00125; 7019.031)
#
#   Source:
#     Lampinen, Jouni, and Zelinka, Ivan (1999).
#     Mechanical engineering design optimization
#     by differential evolution.
#     In: David Corne, Marco Dorigo and Fred Glover (Editors),
#     New Ideas in Optimization, McGraw-Hill, pp 127-146
pressure_vessel <-
    list(obj = function(x) {
             x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]
             0.6224*x1*x3*x4 + 1.7781*x2*x3^2 +
             3.1611*x1^2*x4 + 19.84*x1^2*x3
         },
         con = function(x) {
             x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]
             c(0.0193*x3 - x1,
               0.00954*x3 - x2,
               750.0*1728.0 - pi*x3^2*x4 - 4/3*pi*x3^3)
         })

JDEoptim(c( 1.1,  0.6,   0.0,   0.0),
         c(12.5, 12.5, 240.0, 240.0),
         fn = pressure_vessel$obj, constr = pressure_vessel$con,
         tol = 1e-7, trace = TRUE, triter = 50)
}
