require(DEoptimR)

c.time <- function(...) cat('Time elapsed: ', ..., '\n')
S.time <- function(expr) c.time(system.time(expr))
(doExtras <- DEoptimR:::doExtras() && requireNamespace("mirai", quietly = TRUE))

# Bound-constrained test problems ----------------------------------------------

ap <- function(x) {
#   Aluffi-Pentini's problem
#
#   -10 <= x1, x2 <= 10
#   The function has two local minima, one of them is global
#   with f(x*) ~~ -0.3523 at (-1.0465, 0).
#
#   Source:
#     Ali, M. Montaz, Khompatraporn, Charoenchai, and Zabinsky, Zelda B. (2005).
#     A numerical evaluation of several stochastic algorithms on selected
#     continuous global optimization test problems.
#     Journal of Global Optimization 31, 635-672.

    0.25*x[1]^4 - 0.5*x[1]^2 + 0.1*x[1] + 0.5*x[2]^2
}

set.seed(42)
S.time(ap_seq <- SPJDEoptim(-c(10, 10), c(10, 10), ap,
                            sequential_eval = "objective",
                            tol = 1e-7))

if (doExtras) {
daemons(1, dispatcher = FALSE) # stays within 2-core limit

# slow
set.seed(42)
S.time(ap_par <- SPJDEoptim(-c(10, 10), c(10, 10), ap, tol = 1e-7))
}

# Only inequality constraints --------------------------------------------------

g01 <- list(obj = function(x) {
            #   0 <= xi <= 1 (i = 1, ..., 9), 0 <= xi <= 100 (i = 10, 11, 12)
            #   and 0 <= x13 <= 1
            #   The global minimum is at
            #   x* = (1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1)
            #   where six constraints are active (g1, g2, g3, g7, g8, and g9)
            #   and f(x*) = -15.
            #
            #   Source:
            #     Runarsson, Thomas P., and Yao, Xin (2000).
            #     Stochastic ranking for constrained evolutionary optimization.
            #     IEEE Transactions on Evolutionary Computation 4, 284-294.

                5*(sum(x[1:4]) - crossprod(x[1:4])) - sum(x[5:13])
            },
            con = function(x) {
                x1 <- x[1]; x2 <- x[2]; x3 <- x[3]
                x10 <- x[10]; x11 <- x[11]; x12 <- x[12]
                c(2*(x1 + x2) + x10 + x11 - 10,
                  2*(x1 + x3) + x10 + x12 - 10,
                  2*(x2 + x3) + x11 + x12 - 10,
                  -8*x1 + x10,
                  -8*x2 + x11,
                  -8*x3 + x12,
                  -2*x[4] - x[5] + x10,
                  -2*x[6] - x[7] + x11,
                  -2*x[8] - x[9] + x12)
            })

set.seed(42)
S.time(g01_seq <- SPJDEoptim(rep(0, 13), c(rep(1, 9), rep(100, 3), 1),
                             fn = g01$obj, constr = g01$con,
                             sequential_eval = "both",
                             tol = 1e-7))

if (doExtras) {
# very slow
set.seed(42)
S.time(g01_par <- SPJDEoptim(rep(0, 13), c(rep(1, 9), rep(100, 3), 1),
                             fn = g01$obj, constr = g01$con, tol = 1e-7))

daemons(0) # reset
Sys.sleep(1)
}

# Expected optimal values ------------------------------------------------------
bare_p_v <- function(r) unlist(unname( r[c("par", "value")] ))
stopifnot(
  all.equal( bare_p_v(ap_seq), c(-1.0465, 0, -0.3523),
             tolerance = 1e-3 ),
  all.equal( bare_p_v(g01_seq), c(rep(1, 9), rep(3, 3), 1, -15),
             tolerance = 1e-4 )
)
if (doExtras) {
stopifnot(
  identical( bare_p_v(ap_seq), bare_p_v(ap_par) ),
  identical( bare_p_v(g01_seq), bare_p_v(g01_par) )
)
}

c.time(proc.time())
