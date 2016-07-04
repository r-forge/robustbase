### lmrob()  with "real data"

library(robustbase)

set.seed(0)
data(salinity)
summary(m0.sali  <- lmrob(Y ~ . , data = salinity))
(A1 <- anova(m0.sali, Y ~ X1 + X3))
## -> X2 is not needed
(m1.sali  <- lmrob(Y ~ X1 + X3, data = salinity))
(A2 <- anova(m0.sali, m1.sali)) # the same as before
stopifnot(all.equal(A1[2,"Pr(>chisq)"],
		    A2[2,"Pr(>chisq)"], tolerance=1e-14))
anova(m0.sali, m1.sali, test = "Deviance")
## whereas 'X3' is highly significant:
m2 <- update(m0.sali, ~ . -X3)
(A3 <- anova(m0.sali, m2))
(A4 <- anova(m0.sali, m2, test = "Deviance"))
cX3 <- c(Estimate = -0.627327396, `Std. Error` = 0.15844971,
         `t value` = -3.9591577, `Pr(>|t|)` = 0.000584156)
stopifnot(all.equal(cX3, coef(summary(m0.sali))["X3",], tolerance = 1e-6))


##  example(lmrob)
set.seed(7)
data(coleman)
summary( m1 <- lmrob(Y ~ ., data=coleman) )
stopifnot(c(3,18) == which(m1$w < 0.2))

if(FALSE) # to find out *why setting = "KS201x" fails
trace(lmrob.S, exit = quote({cat("coef:\n"); print(b$coefficients)}))
if(FALSE) # to find out via  setting = "KS201x" fails here in the *initial* estimate
debug(lmrob.S)

data(starsCYG)
  lmST <-    lm(log.light ~ log.Te, data = starsCYG)
(RlmST <- lmrob(log.light ~ log.Te, data = starsCYG, control=lmrob.control(trace = 1)))
summary(RlmST)
## Least Sq. w/ negative slope, where robust has slope ~= 2.2 :
stopifnot(coef(lmST)[["log.Te"]] < 0,
          all.equal(coef(RlmST), c("(Intercept)" = -4.969, log.Te=2.253), tol = 1e-3),
          c(11,20,30,34) == which(RlmST$w < 0.01))
## ==> Now see that  "KS2011" and "KS2014" both break down -- and it is the fault of "lqq" *only* :
(RlmST.11 <- update(RlmST, control = lmrob.control("KS2011",                             trace= 1)))
(RlmST.14 <- update(RlmST, control = lmrob.control("KS2014",                             trace= 1)))
(RlmSTM11 <- update(RlmST, control = lmrob.control("KS2011", method="MM",                trace= 1)))
(RlmSTM14 <- update(RlmST, control = lmrob.control("KS2014", method="MM",                trace= 1)))
## using "biweight" instead of "lqq"  fixes the problem :
(RlmSTM11b <- update(RlmST,control = lmrob.control("KS2011", method="MM", psi="biweight",trace= 1)))
(RlmSTM14b <- update(RlmST,control = lmrob.control("KS2014", method="MM", psi="biweight",trace= 1)))
(RlmST.11b <- update(RlmST,control = lmrob.control("KS2011",              psi="biweight",trace= 1)))
(RlmST.14b <- update(RlmST,control = lmrob.control("KS2014",              psi="biweight",trace= 1)))
## NB: RlmST has component 'init.S' the others have "init" -- documented in ?lmrob.fit == ../man/lmrob.fit.Rd
R.ini.cf <- t(sapply(mget(ls(patt = "^RlmST")), function(OB) OB$init$coef))
R..cf    <- t(sapply(mget(ls(patt = "^RlmST")), coef))
cbind(R.ini.cf, R..cf) ##---> "lqq" is *NOT* robust enough here -- but "biweight" is !!

## Directly look at init.S():
x.s <- model.matrix(~ log.Te, data = starsCYG)
y.s <- model.response(model.frame(log.light ~ log.Te, data = starsCYG))
ini.df <- lmrob.S(x.s, y.s, control=lmrob.control())
ini.11 <- lmrob.S(x.s, y.s, control=lmrob.control("KS2011"))
ini.14 <- lmrob.S(x.s, y.s, control=lmrob.control("KS2014"))
## but these are fine !! :
rbind(deflt = ini.df$coef,
      KS.11 = ini.11$coef,
      KS.14 = ini.14$coef)

##==> No, it is *not* the init.S()


set.seed(47)
data(hbk)
m.hbk <- lmrob(Y ~ ., data = hbk)
summary(m.hbk)
stopifnot(1:10 == which(m.hbk$w < 0.01))

data(heart)
summary(mhrt <- lmrob(clength ~ ., data = heart))
stopifnot(8  == which(mhrt$w < 0.15),
          11 == which(0.61 < mhrt$w & mhrt$w < 0.62),
          c(1:7,9:10,12) == which(mhrt$w > 0.90))

data(stackloss)
mSL <- lmrob(stack.loss ~ ., data = stackloss)
summary(mSL)


cat('Time elapsed: ', proc.time(),'\n') # "stats"
