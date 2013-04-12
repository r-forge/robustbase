library(robustbase)

## From latest  system.file("test-tools-1.R", package="Matrix") :
assert.EQ <- function(target, current, tol = if(show) 0 else 1e-15,
                      show = FALSE, ...) {
    ## Purpose: check equality *and* show non-equality
    ## ----------------------------------------------------------------------
    ## show: if TRUE, return (and hence typically print) all.equal(...)
    if(show) all.equal(target, current, tol = tol)
    else if(!isTRUE(r <- all.equal(target, current, tol = tol)))
	stop("all.equal() |->  ", r)
}

###>> 1 ------------------- family = poisson ------------------------------------

### very simple model [with outliers]
set.seed(113)
y <- rpois(17, lambda = 4)
y[1:2] <- 99:100 # outliers
y
options(digits=10)
rm1 <- glmrob(y ~ 1, family = poisson, trace = TRUE,
              acc = 1e-13) # default is just 1e-4
cm1 <- glm   (y ~ 1, family = poisson, trace = TRUE)

assert.EQ(c(0.0287933850640724, 0.0284930623638766,
		      0.950239140568007, 0.874115394740014),
		    local({w <- rm1$w.r; w[ w != 1 ] }), tol = 1e-14)
assert.EQ(coef(rm1), c("(Intercept)" = 1.41710946076738),tol = 1e-14)

if(FALSE) # for manual digging:
debug(robustbase:::glmrobMqle)

allresid <- function(obj, types = c("deviance", "pearson", "working", "response"))
{
    sapply(types, residuals, object = obj)
}

okFit <- function(obj, check.attr=FALSE, ...) {
  all.equal(obj$y, obj$fitted.values + residuals(obj, "response"),
            check.attr=check.attr, ...)
}

## check validity of several methods simultaneously:
y. <- model.response(model.frame(rm1))
stopifnot(okFit(cm1), okFit(rm1), y. == y)



alr.c <- allresid(cm1)
alr.r <- allresid(rm1)

## MM --- just for now --
plot(resid(cm1),                resid(rm1), asp=1); abline(0,1, col=2)
plot(resid(cm1,type="pearson"), resid(rm1, type="pearson"), asp=1); abline(0,1, col=2)
plot(resid(cm1,type="working"), resid(rm1, type="working"), asp=1); abline(0,1, col=2)

## leave away the outliers --
cm0 <- glm   (y ~ 1, family = poisson, trace = TRUE, subset = -(1:2))
plot(resid(cm0),                resid(rm1)[-(1:2)], asp=1); abline(0,1, col=2)
plot(resid(cm0,type="pearson"), resid(rm1, type="pearson")[-(1:2)], asp=1); abline(0,1, col=2)
plot(resid(cm0,type="working"), resid(rm1, type="working")[-(1:2)], asp=1); abline(0,1, col=2)
plot(resid(cm0,type="response"), resid(rm1, type="response")[-(1:2)], asp=1); abline(0,1, col=2)


## Use weights (slow convergence !)
w2 <- c(rep(1,8), rep(10,9))
rm2 <- glmrob(y ~ 1, family = poisson, trace = TRUE,
              weights = w2, maxit = 500, acc = 1e-10) # default is just 1e-4
## slow convergence
stopifnot(okFit(rm2))


###>> 2 ------------------- family = binomial -----------------------------------

## Using  *factor*  y ...
x <- seq(0,5, length = 120)
summary(px <- plogis(-5 + 2*x))
set.seed(7)
(f <- factor(rbinom(length(x), 1, prob=px)))

summary(m.c0 <- glm   (f ~ x, family = binomial))
summary(m.r0 <- glmrob(f ~ x, family = binomial))

## add outliers --- in y:
f. <- f
f.[i1 <- 2:3] <- 1
f.[i0 <- 110+c(1,7)] <- 0
        m.c1 <- glm   (f. ~ x, family = binomial)
summary(m.r1 <- glmrob(f. ~ x, family = binomial)) ## hmm, not so robust?
stopifnot(m.r1$w.r[c(i0,i1)] < 1/3, # well, at least down weighted
	  ## and coefficients change less :
	  (coef(m.r1) - coef(m.c0)) / (coef(m.c1) - coef(m.c0)) < 1)
assert.EQ(c("(Intercept)" = -3.10817337603974, x = 1.31618564057790),
	  coef(m.r1), tol= 1e-14)

y <- as.numeric(as.character(f.))
m.r2 <- BYlogreg(x0=x, y=y, trace=TRUE)

assert.EQ(c("(Intercept)" = -2.9554950286, x0 = 1.2574679132),
          ## 32-bit{ada-5}  -2.95549502890363   1.25746791332613
	  m.r2$coef, tol=4e-10, show=TRUE)
assert.EQ( c(0.685919891749065, 0.256419206157062),
          ## 32-bit{ada-5}:
          ## 0.685919891858219, 0.256419206203016)
	  m.r2$sterror, tol=4e-10)

if(require("catdata"))
    data(foodstamp) else
load(system.file("external/foodstamp.rda",
		 package="robustbase", mustWork=TRUE))

str(foodstamp) ## y ~ TEN, SUP, INC  {what ugly names!}
m.fsQL <- glmrob(y ~ ., family = binomial, data=foodstamp)
X.fs <- model.matrix(m.fsQL)
stopifnot(dim(X.fs) == c(150, 4)) # including intercept!
m.fsWBY <- BYlogreg(x0=X.fs, y=foodstamp[["y"]],
		    addIntercept=FALSE, trace=TRUE)
m.fs.BY <- BYlogreg(x0=X.fs, y=foodstamp[["y"]], initwml=FALSE,
		    addIntercept=FALSE, trace=TRUE)

