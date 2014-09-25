\name{rlas}
\alias{rlas}
\title{Recursive Robust Median-like Location and Scale}
\description{
  Calculate a recursive estimate of location, asymptotically
  equivalent to the median, and a recursive estimate of scale
  equal to the \bold{MEAN} absolute deviation.
}
\usage{
rlas(y, b = 0.2, mfn = function(n) 0.1 * n^(-0.25),
     nstart = 30, scon = NULL)
}
\arguments{
  \item{y}{numeric vector of i.i.d. data whose location and scale
    parameters are to be estimated.}
  \item{b}{numeric tuning parameter (default value equal to that used by
    Holst, 1987).}
  \item{mfn}{a \code{\link{function}} of the index of the data which
    must be positive and and tend to 0 as the index tends to infinity.
    The default function is that used by Holst, 1987.}
  \item{nstart}{number of starting values: The algorithm is initialized from
    from the first \code{nstart} values of \code{y}.}
  \item{scon}{value for the scale parameter \code{s}.  If provided, the
    scale is held constant throughout and the algorithm becomes
    equivalent to the algorithm of Holst.}
}
\value{
  A list with entries
  \item{locn}{The successive recursive estimates of location.  The
    first \code{nstart - 1} of these are \code{NA}.}
  \item{scale}{The successive recursive estimates of scale. If
    \code{scon} is specified, these are all equal to \code{scon}.}
}
\references{
  Cameron, Murray A. and Turner, T. Rolf (1993).
  Recursive location and scale estimators.
  \emph{Commun. Statist. --- Theory Meth.} \bold{22}(9) 2503--2515.

  Holst, U. (1987).
  Recursive estimators of location.
  \emph{Commun. Statist. --- Theory Meth.} \bold{16} (8) 2201--2226.
}
\author{
 \email{r.turner@auckland.ac.nz}
 \url{http://www.stat.auckland.ac.nz/~rolf}
}
\examples{
set.seed(42)
y  <- rnorm(10000)
z1 <- rlas(y)
z2 <- rlas(y,scon=1)
z3 <- rlas(y,scon=100)
OP <- par(mfrow=c(3,1))
plot(z1$locn,type="l")
abline(h=median(y),col="red")
plot(z2$locn,type="l")
abline(h=median(y),col="red")
plot(z3$locn,type="l")
abline(h=median(y),col="red")
par(OP)
}
\keyword{univar}
\keyword{robust}
