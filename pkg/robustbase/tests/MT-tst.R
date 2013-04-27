require("robustbase")

##---> ./poisson-ex.R
##     ~~~~~~~~~~~~~~  for more glmrobMT() tests

source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
## -> assertError(), showSys.time(), ...
(doExtras <- robustbase:::doExtras())

## The simple intercept example from  ./glmrob-1.R
set.seed(113)
y <- rpois(17, lambda = 4)
y[1:2] <- 99:100 # outliers
y.1 <- y
## To call Victor's version of glmrobMT()  for this case, we need  "hoop jumps"
x.1 <- matrix(0, nrow=length(y.1), ncol=0)

r <- glmrobMT(x.1, y.1, nsubm=100)# some output
str(r)

# was    c(ini = 1.30833281965018, est = 1.29369680430613)
r.64b <- c(ini = 1.30833281965018, est = 1.29369680422799)
stopifnot(r$converged,
	  all.equal(r$initial, r.64b[["ini"]], tol = 1e-13),# rel.diff: 3.394.e-16
	  all.equal(r$final,   r.64b[["est"]], tol = 1e-13) # rel.diff: 2.9178e-15
          )

## now, as the algorithm has a random start:
set.seed(7)
nSim <- if(doExtras) 20 else 2
showSys.time(LL <- replicate(nSim, glmrobMT(x.1, y.1, trace.lev=0),
                             simplify=FALSE))
ini <- sapply(LL, `[[`, "initial")
est <- sapply(LL, `[[`, "final")
## surprise:  all the 20 initial estimators are identical:
stopifnot(diff(range(ini)) == 0,
          diff(range(est)) == 0,
## probably too accurate ... but ok, for now
          all.equal(est[1], r.64b[["est"]], tol = 1e-13),
          all.equal(ini[1], r.64b[["ini"]], tol = 1e-13),
          TRUE)



cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
