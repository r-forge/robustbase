source(system.file("deopt.R", package = "robustbase"))
source("opt-test-funs.R")
set.seed(2345)
sf1. <- jde(c(-100, -100), c(100, 100), 1, sf1, tol = 1e-7)
swf. <- jde(rep(-500, 10), rep(500, 10), 1, swf, tol = 1e-7)
RND. <- jde(c(1e-5, 1e-5), c(16, 16), 1, RND$obj, RND$con, tol = 1e-7)
HEND. <- jde(c(100, 1000, 1000, 10, 10), c(10000, 10000, 10000, 1000, 1000), 1,
             HEND$obj, HEND$con, tol = 0.1)
alkylation. <- jde(c(1500, 1, 3000, 85, 90, 3, 145), c(2000, 120, 3500, 93, 95, 12, 162), 1,
                   alkylation$obj, alkylation$con, tol = 0.1)
stopifnot(
    all.equal( unlist(unname(sf1.[c("par", "value")])), c(0, 0, 0),
               tolerance = 1e-4 ),
    all.equal( unlist(unname(swf.[c("par", "value")])), c(rep(420.97, 10), -418.9829*10),
               tolerance = 1e-5 ),
    all.equal( unlist(unname(RND.[c("par", "value")])), c(3.036504, 5.096052, -0.388812),
               tolerance = 1e-3 ),
    all.equal( unlist(unname(HEND.[c("par", "value")])),
               c(579.19, 1360.13, 5109.92, 182.01, 295.60, 7049.25),
               tolerance = 1e-3 ),
    all.equal( unlist(unname(alkylation.[c("par", "value")])),
               c(1698.256922, 54.274463, 3031.357313, 90.190233,
                 95.0, 10.504119, 153.535355, -1766.36),
               tolerance = 1e-4 ) )