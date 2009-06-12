library(lattice)

if(require("ICS")) {
  # ICS: Compare M-estimate [Max.Lik. of t_{df = 2}] with high-breakdown :
  stopifnot(require("MASS"))

  X <- rbwheel(n, 7, frac = 1/5, spherize=TRUE) # << more than 1/p outliers!
  X.paM <- ics(X, S1=cov, S2= function(.) cov.trob(., nu=2)$cov, stdKurt=FALSE)
  X.paM.<- ics(X, S1=cov, S2= function(.) tM(., df=2)$V, stdKurt = FALSE)
  X.paR <- ics(X, S1=cov, S2= function(.) covMcd(.)$cov, stdKurt = FALSE)

  par.s <- list(plot.symbol = list(alpha=0.4, cex= .5, pch = 16))
  p0 <- splom(~ X, xlab = "X <- rbwheel(500,7, frac = 1/5, spherize=TRUE)",
              pscales = 0,
              par.settings = list(plot.symbol = c(par.s$plot.symbol,
                                  col = "dark gray")))
  p1 <- splom(~ X.paM @ Scores, par.settings=par.s,
              pscales = 0,
              xlab = "ics(X, .. S2= cov.trob(., nu=2)")
  p2 <- splom(~ X.paM.@ Scores, par.settings=par.s,
              pscales = 0,
              xlab = "ics(X, .. S2= tM(., df=2)")
  p3 <- splom(~ X.paR @ Scores, par.settings=par.s,
              pscales = 0,
              xlab = "ics(X, .. S2= covMcd(.)")
  print(p0, split=c(1,1, 2,2), more=TRUE)
  print(p1, split=c(1,2, 2,2), more=TRUE)
  print(p2, split=c(2,1, 2,2), more=TRUE)
  print(p3, split=c(2,2, 2,2))

}

