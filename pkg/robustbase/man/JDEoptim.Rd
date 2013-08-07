\name{JDEoptim}
\alias{JDEoptim}
\title{Constrained Nonlinear Minimization via Differential Evolution}
\usage{
  JDEoptim(lower, upper, f0, fn, constr = NULL, meq = 0,
    slack = 1e-05, NP = 10 * d, Fl = 0.1, Fu = 1,
    tau1 = 0.1, tau2 = 0.1, tau3 = 0.1,
    jitter_factor = 0.001, tol = 1e-15,
    criteria = c("median", "max"), maxiter = 200 * d,
    init = NULL, trace = FALSE, triter = 1,
    verbose = FALSE, ...)
}
\arguments{
  \item{lower}{}
  \item{upper}{}

  \item{f0}{}

  \item{fn}{the \code{\link{function}} to be minimized.}

  \item{constr}{\code{NULL} or \code{\link{function}} for specifying
  (nonlinear) constraints.}

  \item{meq}{}

  \item{slack}{}

  \item{NP}{}

  \item{Fl}{}

  \item{Fu}{}

  \item{tau1}{}

  \item{tau2}{}

  \item{tau3}{}

  \item{jitter_factor}{}

  \item{tol}{}

  \item{criteria}{function, or character string specifying
  a function which is applied to the population and
  compared with the so-far best one ....FIXME}

  \item{maxiter}{}

  \item{init}{}

  \item{trace}{}

  \item{triter}{}

  \item{verbose}{}

  \item{...}{optionially arguments passed to \code{fn()}
  \emph{and} \code{constr()} if that is not \code{NULL}}
}
\description{
  Constrained nonlinear minimization (Differential
  Evolution [1]). It uses the DE/rand/1/either-or mutation
  strategy [2]. It also utilizes only one population [3]
  against two sets as in the original DE algorithm [4].
}
\details{
  Copyright 2013, Eduardo L. T. Conceicao Available under
  the GPL (>= 2)
}
\value{
  an object ... FIXME
}
\author{
  Eduardo L. T. Conceicao
}
\references{
  [1] Brest, J., Greiner, S., Boskovic, B., Mernik, M. and Zumer, V. (2006)
  Self-adapting control parameters in differential evolution: a
  comparative study on numerical benchmark problems.
  \emph{IEEE Trans. Evol. Comput.} \bold{10}, 646--657.

  [2] Price, Kenneth V., Storn, Rainer M., and Lampinen, Jouni A. (2005)
  \emph{Differential Evolution: A practical approach to global optimization}.
  Springer, Berlin, pp. 117--118.

  [3] Babu, B. V., and Angira, Rakesh (2006)
  Modified differential evolution (MDE) for optimization of
  non-linear chemical processes.
  \emph{Computers and Chemical Engineering} \bold{30}, 989--1002.

  [4] Storn, Rainer, and Price, Kenneth (1997)
  Differential evolution - a simple and efficient heuristic for global
  optimization over continuous spaces.
  \emph{Journal of Global Optimization} \bold{11}, 341--359.
}
\seealso{
  Function \code{\link[DEoptim]{DEoptim}} in the
  \pkg{DEoptim} package has many more options than
  \code{JDEoptim()}, but does not allow constraints in the same
  flexible manner.
}
\examples{
source(system.file("xtraR/opt-test-funs.R", package = "robustbase"))
## --> sf1(), swf() functions +
##  RND, HEND, and alkylation, each list of $obj and $con  testing functions

sf1 # the function to minimize
JDEoptim(c(-100, -100), c(100, 100),  1, sf1, tol = 1e-7, trace = TRUE)

swf # potentially high-dim. function to minimize
JDEoptim(rep(-500, 10), rep(500, 10), 1, swf, tol = 1e-7, trace = TRUE)

JDEoptim(c(1e-5, 1e-5), c(16, 16), 1,
    RND$obj, RND$con, tol = 1e-7, trace = TRUE)
JDEoptim(c(100, 1000, 1000, 10, 10), c(10000, 10000, 10000, 1000, 1000), 1,
    HEND$obj, HEND$con, tol = 1e-7, trace = TRUE)
JDEoptim(c(1500, 1, 3000, 85, 90, 3, 145), c(2000, 120, 3500, 93, 95, 12, 162), 1,
    alkylation$obj, alkylation$con, tol = 1e-7, trace = TRUE)
}
