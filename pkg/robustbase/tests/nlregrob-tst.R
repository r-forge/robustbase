stopifnot(require("robustbase"))
## testing functions:
source(system.file("xtraR/ex-funs.R", package = "robustbase"))
c.time <- function(...) cat('Time elapsed: ', ..., '\n')
S.time <- function(expr) c.time(system.time(expr))
showProc.time <- local({ ## function + 'pct' variable
    pct <- proc.time()
    function(final="\n") { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	## 'Time ..' *not* to be translated:  tools::Rdiff() skips its lines!
	cat('Time elapsed: ', (pct - ot)[1:3], final)
    }
})

## as long as we don't export these (nor provide an nlrob(., method=.) interface:
nlrob.MM  <- robustbase:::nlrob.MM
nlrob.tau <- robustbase:::nlrob.tau
nlrob.CM  <- robustbase:::nlrob.CM
nlrob.mtl <- robustbase:::nlrob.mtl

(doExtras <- robustbase:::doExtras())
if(doExtras) {
    NP <- 20
} else { ## fast
    NP <- 10
}

use.true <- !doExtras # (but not necessarily ..)
if(use.true) { # population size = NP (random) + 1 (true parameters)
    true_params       <- c(1, 0.2)
    true_params_sigma <- c(1, 0.2, 1)
} else {
    true_params <- true_params_sigma <- NULL
}

## Stromberg, Arnold J. (1993).
## Computation of high breakdown nonlinear regression parameters.
## J. Amer. Statist. Assoc. 88(421), 237-244.

## exponential regression
Expo <- function(x, a, b) exp(a + b*x)
set.seed(2345) # for reproducibility
## data without outliers:
d.exp30 <- data.frame(x = sort( runif(30, 0, 10) ), err = rnorm(30))
d.exp30 <- transform(d.exp30, y = Expo(x, 1, 0.2) + err)
## classical
Cfit <- nls(y ~ Expo(x, a, b), data = d.exp30, start = c(a = 1, b = 0.2),
            control = nls.control(tol = 5e-8, printEval = TRUE))
showProc.time()

## robust
Rfit.MM.S.bisquare <-
    nlrob.MM( y ~ Expo(x, a, b), data = d.exp30, pnames = c("a", "b"),
              lower = c(-10, -2), upper = c(10, 2),
              NP = NP, add_to_init_pop = true_params )
Rfit.MM.S.lqq        <- update(Rfit.MM.S.bisquare, psi = "lqq")
Rfit.MM.S.optimal    <- update(Rfit.MM.S.bisquare, psi = "optimal")
Rfit.MM.S.hampel     <- update(Rfit.MM.S.bisquare, psi = "hampel")
showProc.time()
Rfit.MM.lts.bisquare <- update(Rfit.MM.S.bisquare, init = "lts")
Rfit.MM.lts.lqq      <- update(Rfit.MM.S.bisquare, init = "lts", psi = "lqq")
Rfit.MM.lts.optimal  <- update(Rfit.MM.S.bisquare, init = "lts", psi = "optimal")
Rfit.MM.lts.hampel   <- update(Rfit.MM.S.bisquare, init = "lts", psi = "hampel")
showProc.time()

S.time(Rfit.tau.bisquare <-
    nlrob.tau( y ~ Expo(x, a, b), data = d.exp30, pnames = c("a", "b"),
               lower = c(-10, -2), upper = c(10, 2),
               NP = NP, add_to_init_pop = true_params ))
S.time(Rfit.tau.optimal <- update(Rfit.tau.bisquare, psi = "optimal"))

S.time(Rfit.CM <- nlrob.CM( y ~ Expo(x, a, b), data = d.exp30,
                            pnames = c("a", "b", "sigma"),
                            lower = c(-10, -2, 0), upper = c(10, 2, 10),
                            NP = NP, add_to_init_pop = true_params_sigma ))
S.time(Rfit.mtl <- nlrob.mtl( y ~ Expo(x, a, b), data = d.exp30,
                              pnames = c("a", "b", "sigma"),
                              lower = c(-10, -2, 0), upper = c(10, 2, 3),
                              NP = NP, add_to_init_pop = true_params_sigma ))
showProc.time()

## 40% outliers present {use different data name: seen in print(<fitted model>)
d.exp40out <- within(d.exp30, y[15:27] <- y[15:27] + 100)

Rf.out.MM.S.bisquare   <- update(Rfit.MM.S.bisquare, data=d.exp40out)
Rf.out.MM.S.lqq        <- update(Rf.out.MM.S.bisquare, psi = "lqq")
Rf.out.MM.S.optimal    <- update(Rf.out.MM.S.bisquare, psi = "optimal")
Rf.out.MM.S.hampel     <- update(Rf.out.MM.S.bisquare, psi = "hampel")
showProc.time()
Rf.out.MM.lts.bisquare <- update(Rf.out.MM.S.bisquare, init= "lts")
Rf.out.MM.lts.lqq      <- update(Rf.out.MM.S.bisquare, init= "lts", psi= "lqq")
Rf.out.MM.lts.optimal  <- update(Rf.out.MM.S.bisquare, init= "lts", psi="optimal")
Rf.out.MM.lts.hampel   <- update(Rf.out.MM.S.bisquare, init= "lts", psi= "hampel")
showProc.time()

Rf.out.tau.bisquare <- update(Rfit.tau.bisquare, data=d.exp40out)
Rf.out.tau.optimal  <- update(Rfit.tau.bisquare, data=d.exp40out, psi = "optimal")

Rf.out.CM  <- update(Rfit.CM,  data=d.exp40out)
Rf.out.mtl <- update(Rfit.mtl, data=d.exp40out)
showProc.time()

## presence of high leverage point outliers
d.exp.Hlev <- within(d.exp40out, {
    x[28:30] <- x[28:30] + 10   ## shift  10
    y <- Expo(x, 1, 0.2) + err
    y[28:30] <- y[28:30] + 500
})

Rf.Hlev.MM.S.bisquare   <- update(Rfit.MM.S.bisquare, data = d.exp.Hlev)
Rf.Hlev.MM.S.lqq        <- update(Rf.Hlev.MM.S.bisquare, psi = "lqq")
Rf.Hlev.MM.S.optimal    <- update(Rf.Hlev.MM.S.bisquare, psi = "optimal")
Rf.Hlev.MM.S.hampel     <- update(Rf.Hlev.MM.S.bisquare, psi = "hampel")
showProc.time()
Rf.Hlev.MM.lts.bisquare <- update(Rf.Hlev.MM.S.bisquare, init="lts")
Rf.Hlev.MM.lts.lqq      <- update(Rf.Hlev.MM.S.bisquare, init="lts", psi= "lqq")
Rf.Hlev.MM.lts.optimal  <- update(Rf.Hlev.MM.S.bisquare, init="lts", psi="optimal")
Rf.Hlev.MM.lts.hampel   <- update(Rf.Hlev.MM.S.bisquare, init="lts", psi= "hampel")
showProc.time()

Rf.Hlev.tau.bisquare <- update(Rfit.tau.bisquare, data = d.exp.Hlev)
Rf.Hlev.tau.optimal  <- update(Rf.Hlev.tau.bisquare, psi = "optimal")

Rf.Hlev.CM  <- update(Rfit.CM,  data = d.exp.Hlev)
Rf.Hlev.mtl <- update(Rfit.mtl, data = d.exp.Hlev)
showProc.time()

				        cfcl <- coef(Cfit)
## no outliers present
assert.EQ(coef(Rfit.MM.S.bisquare),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.S.lqq),		cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.S.optimal),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.S.hampel),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.lts.bisquare),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.lts.lqq),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.lts.optimal),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.MM.lts.hampel),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.tau.bisquare),	cfcl, tol = 0.02, giveRE=TRUE)# 0.009873
assert.EQ(coef(Rfit.tau.optimal),	cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.CM)[-3],		cfcl, tol = 0.01, giveRE=TRUE)
assert.EQ(coef(Rfit.mtl)[-3],		cfcl, tol = 0.05, giveRE=TRUE)
## 40% outliers present
assert.EQ(coef(Rf.out.MM.S.bisquare),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.S.lqq),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.S.optimal),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.S.hampel),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.lts.bisquare),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.lts.lqq),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.lts.optimal),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.MM.lts.hampel),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.tau.bisquare),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.tau.optimal),	cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.CM)[-3],		cfcl, tol = 0.1, giveRE=TRUE)
assert.EQ(coef(Rf.out.mtl)[-3],		cfcl, tol = 0.1, giveRE=TRUE)
## presence of high leverage point outliers
assert.EQ(coef(Rf.Hlev.MM.S.bisquare),	cfcl, tol = .1, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.S.lqq),	cfcl, tol = .1, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.S.optimal),	cfcl, tol = .1, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.S.hampel),	cfcl, tol = .1, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.lts.bisquare),cfcl,tol = .1, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.lts.lqq),	cfcl, tol = .1, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.lts.optimal), cfcl, tol = .1, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.MM.lts.hampel),	cfcl, tol = .1, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.tau.bisquare),	cfcl, tol = .2, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.tau.optimal),	cfcl, tol = .2, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.CM)[-3],		cfcl, tol = .2, giveRE=TRUE)
assert.EQ(coef(Rf.Hlev.mtl)[-3],	cfcl, tol = .2, giveRE=TRUE)

mods <- sapply(ls.str(patt="^Rf"), get, simplify=FALSE)
prblm <- mods[!sapply(mods, `[[`, "status") != "converged"]
if(length(prblm)) {
    cat("\n*** NON-converged model fits:\n")
    print(prblm)
}
