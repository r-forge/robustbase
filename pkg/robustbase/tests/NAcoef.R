## test handing of NA coefficients / singular fits
## also check:
## -- what would have to be done if class "lm" was added.
## -- general compatibility to class lm.
require(robustbase)

## generate simple example data
data <- expand.grid(x1=letters[1:3], x2=LETTERS[1:3], rep=1:3)
set.seed(1)
data$y <- rnorm(nrow(data))
## drop all combinations of one interation:
data <- subset(data, x1 != 'c' | (x2 != 'B' & x2 != 'C'))
## add collinear variables
data$x3 <- rnorm(nrow(data))
data$x4 <- rnorm(nrow(data))
data$x5 <- data$x3 + data$x4
## add some NA terms
data$y[1] <- NA
data$x4[2:3] <- NA ## to test anova

m0 <- lm(y ~ x1*x2 + x3, data)
m1 <- lm(y ~ x1*x2 + x3 + x4 + x5, data)
set.seed(2)
m2 <- lmrob(y ~ x1*x2 + x3 + x4 + x5, data)
m3 <- lmrob(y ~ x1*x2 + x3 + x4, data)
m4 <- lmrob(y ~ x1*x2 + x3, data)

## clean version of m2 (to check predict)
data2 <- data.frame(y=data$y[-(1:3)], m2$x[,!is.na(m2$coef)])
set.seed(2)
m2c <- lmrob(y ~ x1b + x1c + x2B + x2C + x3 + x4 + x1b:x2B + x1b:x2C , data2)

## add class lm to m2 (for now)
class(m2) <- c(class(m2), "lm")
class(m4) <- c(class(m4), "lm")

## the full matrix (data) should be returned by model matrix (frame)
stopifnot(all.equal(model.matrix(m1), model.matrix(m2)),
          all.equal(model.frame(m1), model.frame(m2)))
## qr decomposition should be for the full data and pivots identical lm result
qr.m1 <- qr(m1)$qr
qr.m2 <- m2$qr$qr
stopifnot(NCOL(qr.m2) == NCOL(qr.m1),
          NROW(qr.m2) == NROW(qr.m1),
          length(m2$qr$qraux) == length(qr(m1)$qraux),
          all.equal(m2$qr$pivot, qr(m1)$pivot),
          all.equal(dimnames(qr.m2),dimnames(qr.m1)))
## the alias function should return the same result
stopifnot(all.equal(alias(m1), alias(m2)))

####
## these helper functions should print NAs for the dropped coefficients
print(m2)
summary(m2)
confint(m2)
## drop1 should return df = 0
#drop1(m2) ## drop.lm does not return valid results (yet)!

####
## methods that should just drop the NA coefficients
## m3 is actually the same as m2, so anova should raise an error
stopifnot(is(suppressWarnings(try(anova(m2, m3, test="Wald"),silent=TRUE)), "try-error"))
stopifnot(is(suppressWarnings(try(anova(m2, m3, test="Deviance"),silent=TRUE)), "try-error"))
## but comparing m2 and m4 should be ok
anova(m2, m4, test="Wald")
anova(m2, m4, test="Deviance")
## commands with single #:
## they do (or might) not return sensible results for robust fits
## and need to be checked again
#cooks.distance(m2)
#deviance(m2)
#dfbeta(m2)
#dfbetas(m2)
#effects(m2) ## fails
#extractAIC(m2)
#stopifnot(all.equal(hatvalues(m2), robustbase:::lmrob.leverages(wqr=m2$qr))) ## fails
#influence(m2)
kappa(m2)
stopifnot(all.equal(labels(m2), labels(m1)))
#logLik(m2)
## plot(m2, which=1) ## fails
par(mfrow=c(2,2))
plot(m2, which=2:4)
stopifnot(all.equal(predict(m2), predict(m2c)))
stopifnot(all.equal(predict(m2, se.fit=TRUE, interval="confidence"),
                    predict(m2c, se.fit=TRUE, interval="confidence")))
predict(m2, type="terms", se.fit=TRUE, interval="confidence")
#proj(m2) ## fails
residuals(m2)
#rstandard(m2)
#rstudent(m2)
#simulate(m2) ## just $weights needs to be changed to prior weights
stopifnot(all.equal(variable.names(m2), variable.names(m1)))
vcov(m2)
stopifnot(all.equal(dimnames(m2), dimnames(m1)))
stopifnot(all.equal(dummy.coef(m2), dummy.coef(m1), tol=1)) ## tol=1 to check only structure

## other helper functions
stopifnot(all.equal(case.names(m2), case.names(m1)),
          all.equal(family(m2), family(m1)),
          all.equal(formula(m2), formula(m1)),
          all.equal(nobs(m2), nobs(m1)))
#add1(m4, ~ . + x3 + x4 + x5) ## does not return valid results (yet)!


## test other initial estimators
lmrob(y ~ x1*x2 + x3 + x4 + x5, data, init="M-S")
lmrob(y ~ x1*x2 + x3 + x4 + x5, data, init=lmrob.lar)
