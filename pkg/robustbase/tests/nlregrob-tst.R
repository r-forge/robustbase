stopifnot(require("robustbase"))
## as long as we don't export these (nor provide an nlrob(., method=.) interface:
nlrob.MM  <- robustbase:::nlrob.MM
nlrob.tau <- robustbase:::nlrob.tau
nlrob.CM  <- robustbase:::nlrob.CM
nlrob.mtl <- robustbase:::nlrob.mtl

fast <- !robustbase:::doExtras()
true_params <- if (fast) c(1, 0.2) else NULL
true_params_sigma <- if (fast) c(1, 0.2, 1) else NULL
NP <- if (fast) 9 else 20 # population size = 9 (random) + 1 (true parameters)
## presence of high leverage point outliers:
## MM-estimates with initial S-estimator converge when the value of the
## true parameters is added to the initial population
excl <- function(x) if (fast) x else TRUE

## Stromberg, Arnold J. (1993).
## Computation of high breakdown nonlinear regression parameters.
## J. Amer. Statist. Assoc. 88(421), 237-244.

set.seed(2345) # for reproducibility
df <- data.frame(x = sort( runif(30, 0, 10) ), err = rnorm(30))
## exponential regression
Expo <- function(x, a, b) exp(a + b*x)

## no outliers present
df <- transform(df, y = Expo(x, 1, 0.2) + err)
## classical
Cfit <- nls(y ~ Expo(x, a, b), data = df, start = c(a = 1, b = 0.2),
            control = nls.control(tol = 5e-8, printEval = TRUE))

## robust
Rfit.MM.S.bisquare <-
    nlrob.MM( y ~ Expo(x, a, b), data = df, pnames = c("a", "b"),
              lower = c(-10, -2), upper = c(10, 2),
              NP = NP, add_to_init_pop = true_params )
Rfit.MM.S.lqq <- update(Rfit.MM.S.bisquare, psi = "lqq")
Rfit.MM.S.optimal <- update(Rfit.MM.S.bisquare, psi = "optimal")
Rfit.MM.S.hampel <- update(Rfit.MM.S.bisquare, psi = "hampel")
Rfit.MM.lts.bisquare <- update(Rfit.MM.S.bisquare, estim = "lts")
Rfit.MM.lts.lqq <- update(Rfit.MM.S.bisquare, estim = "lts", psi = "lqq")
Rfit.MM.lts.optimal <- update(Rfit.MM.S.bisquare, psi = "optimal")
Rfit.MM.lts.hampel <- update(Rfit.MM.S.bisquare, estim = "lts", psi = "hampel")

Rfit.tau.bisquare <-
    nlrob.tau( y ~ Expo(x, a, b), data = df, pnames = c("a", "b"),
               lower = c(-10, -2), upper = c(10, 2),
               NP = NP, add_to_init_pop = true_params )
Rfit.tau.optimal <- update(Rfit.tau.bisquare, psi = "optimal")

Rfit.CM <- nlrob.CM( y ~ Expo(x, a, b), data = df,
                     pnames = c("a", "b", "sigma"),
                     lower = c(-10, -2, 0), upper = c(10, 2, 10),
                     NP = NP, add_to_init_pop = true_params_sigma )
Rfit.mtl <- nlrob.mtl( y ~ Expo(x, a, b), data = df,
                       pnames = c("a", "b", "sigma"),
                       lower = c(-10, -2, 0), upper = c(10, 2, 3),
                       NP = NP, add_to_init_pop = true_params_sigma )

## 40% outliers present
df$y[15:27] <- df$y[15:27] + 100
Rfit.out.MM.S.bisquare <- update(Rfit.MM.S.bisquare)
Rfit.out.MM.S.lqq <- update(Rfit.MM.S.bisquare, psi = "lqq")
Rfit.out.MM.S.optimal <- update(Rfit.MM.S.bisquare, psi = "optimal")
Rfit.out.MM.S.hampel <- update(Rfit.MM.S.bisquare, psi = "hampel")
Rfit.out.MM.lts.bisquare <- update(Rfit.MM.S.bisquare, estim = "lts")
Rfit.out.MM.lts.lqq <- update(Rfit.MM.S.bisquare, estim = "lts", psi = "lqq")
Rfit.out.MM.lts.optimal <-
    update(Rfit.MM.S.bisquare, estim = "lts", psi = "optimal")
Rfit.out.MM.lts.hampel <-
    update(Rfit.MM.S.bisquare, estim = "lts", psi = "hampel")

Rfit.out.tau.bisquare <- update(Rfit.tau.bisquare)
Rfit.out.tau.optimal <- update(Rfit.tau.bisquare, psi = "optimal")

Rfit.out.CM <- update(Rfit.CM)
Rfit.out.mtl <- update(Rfit.mtl)

## presence of high leverage point outliers
## shift = 10
df$x[28:30] <- df$x[28:30] + 10
df <- transform(df, y = Expo(x, 1, 0.2) + err)
df$y[28:30] <- df$y[28:30] + 500
Rfit.leverage.MM.S.bisquare <- update(Rfit.MM.S.bisquare)
Rfit.leverage.MM.S.lqq <- update(Rfit.MM.S.bisquare, psi = "lqq")
Rfit.leverage.MM.S.optimal <- update(Rfit.MM.S.bisquare, psi = "optimal")
Rfit.leverage.MM.S.hampel <- update(Rfit.MM.S.bisquare, psi = "hampel")
Rfit.leverage.MM.lts.bisquare <- update(Rfit.MM.S.bisquare, estim = "lts")
Rfit.leverage.MM.lts.lqq <-
    update(Rfit.MM.S.bisquare, estim = "lts", psi = "lqq")
Rfit.leverage.MM.lts.optimal <-
    update(Rfit.MM.S.bisquare, estim = "lts", psi = "optimal")
Rfit.leverage.MM.lts.hampel <-
    update(Rfit.MM.S.bisquare, estim = "lts", psi = "hampel")

Rfit.leverage.tau.bisquare <- update(Rfit.tau.bisquare)
Rfit.leverage.tau.optimal <- update(Rfit.tau.bisquare, psi = "optimal")

Rfit.leverage.CM <- update(Rfit.CM)
Rfit.leverage.mtl <- update(Rfit.mtl)

stopifnot(
    ## no outliers present
    all.equal( coef(Rfit.MM.S.bisquare), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.MM.S.lqq), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.MM.S.optimal), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.MM.S.hampel), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.MM.lts.bisquare), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.MM.lts.lqq), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.MM.lts.optimal), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.MM.lts.hampel), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.tau.bisquare), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.tau.optimal), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.CM)[-3], coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.mtl)[-3], coef(Cfit), tolerance = 0.1 ),
    ## 40% outliers present
    all.equal( coef(Rfit.out.MM.S.bisquare), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.MM.S.lqq), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.MM.S.optimal), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.MM.S.hampel), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.MM.lts.bisquare), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.MM.lts.lqq), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.MM.lts.optimal), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.MM.lts.hampel), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.tau.bisquare), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.tau.optimal), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.CM)[-3], coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.out.mtl)[-3], coef(Cfit), tolerance = 0.1 ),
    ## presence of high leverage point outliers
    excl( all.equal( coef(Rfit.leverage.MM.S.bisquare), coef(Cfit), tol = 0.1 ) ),
    excl( all.equal( coef(Rfit.leverage.MM.S.lqq), coef(Cfit), tol = 0.1 ) ),
    excl( all.equal( coef(Rfit.leverage.MM.S.optimal), coef(Cfit), tol = 0.1 ) ),
    excl( all.equal( coef(Rfit.leverage.MM.S.hampel), coef(Cfit), tol = 0.1 ) ),
    all.equal( coef(Rfit.leverage.MM.lts.bisquare), coef(Cfit), tol = 0.1 ),
    all.equal( coef(Rfit.leverage.MM.lts.lqq), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.leverage.MM.lts.optimal), coef(Cfit), tol = 0.1 ),
    all.equal( coef(Rfit.leverage.MM.lts.hampel), coef(Cfit), tolerance = 0.1 ),
    all.equal( coef(Rfit.leverage.tau.bisquare), coef(Cfit), tolerance = 0.2 ),
    all.equal( coef(Rfit.leverage.tau.optimal), coef(Cfit), tolerance = 0.2 ),
    all.equal( coef(Rfit.leverage.CM)[-3], coef(Cfit), tolerance = 0.2 ),
    all.equal( coef(Rfit.leverage.mtl)[-3], coef(Cfit), tolerance = 0.5 ) )

