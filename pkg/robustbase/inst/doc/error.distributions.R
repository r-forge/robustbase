## error.distributions.R: additional error distributions for use in simulations

require(skewt)

####################################################################################
## contaminated normal
####################################################################################

rcnorm <- function (n,mean=0,sd=1,epsilon=0.1,meanc=mean,sdc=sqrt(10)*sd) {
  e <- rnorm(n,mean,sd)
  nc <- floor(epsilon*n)
  idx <- sample(1:n,nc)
  e[idx] <- rnorm(nc,meanc,sdc)
  e
}

## ignore other arguments for the moment
pcnorm <- function(q,mean=0,sd=1,lower.tail=TRUE,log.p=FALSE,...)
  pnorm(q,mean,sd,lower.tail,log.p)

## ignore other arguments for the moment
qcnorm <- function(p,mean=0,sd=1,lower.tail=TRUE,log.p=FALSE,...)
  qnorm(p,mean,sd,lower.tail,log.p)

## ignore other arguments for the moment
dcnorm <- function(x,mean=0,sd=1,log=FALSE,...)
  dnorm(x,mean,sd,log)
