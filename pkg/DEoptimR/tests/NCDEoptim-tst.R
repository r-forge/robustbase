require(DEoptimR)

c.time <- function(...) cat('Time elapsed: ', ..., '\n')
S.time <- function(expr) c.time(system.time(expr))
(doExtras <- DEoptimR:::doExtras())

set.seed(2345)
# Bound-constrained test problems ----------------------------------------------

bl <- function(x) {
#   Becker and Lago problem
#
#   -10 <= x1, x2 <= 10
#   The function has four minima located at (+-5, +-5), all with f(x*) = 0.
#
#   Source:
#     Ali, M. Montaz, Khompatraporn, Charoenchai, and Zabinsky, Zelda B. (2005).
#     A numerical evaluation of several stochastic algorithms on selected
#     continuous global optimization test problems.
#     Journal of Global Optimization 31, 635-672.

    sum((abs(x) - 5)^2)
}

S.time(bl_ <- NCDEoptim(-c(10, 10), c(10, 10), bl,
                        niche_radius = 5, maxiter = 100))

# Expected optimal values ------------------------------------------------------

stopifnot(
  all.equal( as.vector(abs(bl_$solution_arch)), rep(5, 8), tolerance = 1e-4 )
)

c.time(proc.time())
