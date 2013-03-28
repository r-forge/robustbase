\name{rbwheel}
\alias{rbwheel}
\title{Multivariate Barrow Wheel Distribution Random Vectors}
\description{
  Generate \eqn{p}-dimensional random vectors according to Stahel's
  Barrow Wheel Distribution.
}
\usage{
rbwheel(n, p, frac = 1/p, sig1 = 0.05, sig2 = 1/10,
        rGood = rnorm,
        rOut = function(n) sqrt(rchisq(n, p - 1)) * sign(runif(n, -1, 1)),
        U1 = rep(1, p),
        scaleAfter = TRUE, scaleBefore = FALSE, spherize = FALSE,
        fullResult = FALSE)
}
\arguments{
  \item{n}{integer, specifying the sample size.}
  \item{p}{integer, specifying the dimension (aka number of variables).}
  \item{frac}{numeric, the proportion of outliers.
    The default, \eqn{1/p}, corresponds to the (asymptotic) breakdown
    point of M-estimators.}
  \item{sig1}{thickness of the \dQuote{wheel}, (\eqn{= \sigma}
    \code{(good[,1])}), a non-negative numeric.}
  \item{sig2}{thickness of the \dQuote{axis} (compared to 1).}
  \item{rGood}{function; the generator for \dQuote{good} observations.}
  \item{rOut}{function, generating the outlier observations.}
  \item{U1}{p-vector to which \eqn{(1,0,\dots,0)} is rotated.}
  \item{scaleAfter}{logical indicating if the matrix is re-scaled \emph{after}
    rotation (via \code{\link{scale}()})..  Default \code{TRUE}; note
    that this used to be false by default in the first public version.}
  \item{scaleBefore}{logical indicating if the matrix is re-scaled before
    rotation (via \code{\link{scale}()}).}
  \item{spherize}{logical indicating if the matrix is to be
    \dQuote{spherized}, i.e., rotated and scaled to have empirical
    covariance \eqn{I_p}.  This means that the principal components are
    used (before rotation).}
  \item{fullResult}{logical indicating if in addition to the \eqn{n
      \times p}{n x p} matrix, some intermediate quantities are returned as
    well.}
}
\details{
 ....  %% FIXME, maybe only  "see the reference."
}
\value{
  By default (when \code{fullResult} is \code{FALSE}), an
  \eqn{n \times p}{n x p} matrix of \eqn{n} sample vectors of the
  \eqn{p} dimensional barrow wheel distribution, with an attribute,
  \code{n1} specifying the exact number of \dQuote{good} observations,
  \eqn{n1 \approx (1-f)\cdot n}{n1 ~= (1-f)*n},  \eqn{f = }\code{frac}.

  If \code{fullResult} is \code{TRUE}, a list with components
  \item{X}{the \eqn{n \times p}{n x p} matrix of above,
    \code{X = X0 \%*\% A}, where \code{A <- \link{Qrot}(p, u = U1)}, and
    \code{X0} is the corresponding matrix before rotation, see below.}
  \item{X0}{.........}
  \item{A}{the \eqn{p \times p}{p x p} rotation matrix, see above.}
  \item{n1}{the number of \dQuote{good} observations, see above.}
  \item{n2}{the number of \dQuote{outlying} observations, \eqn{n2 = n - n1}.}
}
\references{
  \url{http://stat.ethz.ch/people/maechler/robustness}
  %% now contains the talk @ ICORS 2009 (Parma)

  Stahel, W.~A. and \enc{Mächler}{Maechler}, M. (2009).
  Comment on ``invariant co-ordinate selection'',
  \emph{Journal of the Royal Statistical Society B} \bold{71}, 584--586.
}
\author{Werner Stahel and Martin Maechler}
% \seealso{ %............
% }
\examples{
set.seed(17)
rX8 <- rbwheel(1000,8, fullResult = TRUE, scaleAfter=FALSE)
with(rX8, stopifnot(all.equal(X, X0 \%*\% A,    tol = 1e-15),
                    all.equal(X0, X \%*\% t(A), tol = 1e-15)))
##--> here, don't need to keep X0 (nor A, since that is Qrot(p))

## for n = 100,  you  don't see "it", but may guess .. :
n <- 100
pairs(r <- rbwheel(n,6))
n1 <- attr(r,"n1") ; pairs(r, col=1+((1:n) > n1))

## for n = 500, you *do* see it :
n <- 500
pairs(r <- rbwheel(n,6))
## show explicitly
n1 <- attr(r,"n1") ; pairs(r, col=1+((1:n) > n1))

## but increasing sig2 does help:
pairs(r <- rbwheel(n,6, sig2 = .2))

## show explicitly
n1 <- attr(r,"n1") ; pairs(r, col=1+((1:n) > n1))

set.seed(12)
pairs(X <- rbwheel(n, 7, spherize=TRUE))
colSums(X) # already centered

if(require("ICS")) {
  # ICS: Compare M-estimate [Max.Lik. of t_{df = 2}] with high-breakdown :
  stopifnot(require("MASS"))
  X.paM <- ics(X, S1 = cov, S2 = function(.) cov.trob(., nu=2)$cov, stdKurt = FALSE)
  X.paM.<- ics(X, S1 = cov, S2 = function(.) tM(., df=2)$V, stdKurt = FALSE)
  X.paR <- ics(X, S1 = cov, S2 = function(.) covMcd(.)$cov, stdKurt = FALSE)
  plot(X.paM) # not at all clear
  plot(X.paM.)# ditto
  plot(X.paR)# very clear
}
## Similar such experiments --->  demo(rbwheel_d)  and   demo(rbwheel_ics)
##                                --------------         -----------------
}
\keyword{distribution}
\keyword{robust}
