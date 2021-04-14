\name{mvBACON}
\Rdversion{1.1}
\alias{mvBACON}
\title{BACON: Blocked Adaptive Computationally-Efficient Outlier Nominators}
\description{
  This function performs an outlier identification
  algorithm to the data in the x array [n x p] and y vector [n]
  following the lines described by Hadi et al. for their
  BACON outlier procedure.
}
\usage{
mvBACON(x, collect = 4, m = min(collect * p, n * 0.5), alpha = 0.05,
        init.sel = c("Mahalanobis", "dUniMedian", "random", "manual", "V2"),
        man.sel, maxsteps = 100, allowSingular = FALSE, verbose = TRUE)
}
\arguments{
  \item{x}{numeric matrix (of dimension \eqn{[n x p]}),
    not supposed to contain missing values.}
  \item{collect}{a multiplication factor \eqn{c}, when \code{init.sel} is not
    \code{"manual"}, to define \eqn{m}, the size of the initial basic
    subset, as \eqn{m := c \cdot p}{c * p}, in practice,
    \code{m <- min(p * collect, n/2)}.}
  \item{m}{integer in \code{1:n} specifying the \emph{size} of the initial basic
    subset; used only when \code{init.sel} is not \code{"manual"}.}
  \item{alpha}{determines the cutoff value for the Mahalanobis
    distances (see details).}
  \item{init.sel}{character string, specifying the initial selection
    mode; implemented modes are:
    \describe{
      \item{"Mahalanobis"}{based on Mahalanobis distances (default); the
	version \eqn{V1} of the reference; affine invariant but not robust.}
      \item{"dUniMedian"}{based on the distances from the
	\bold{uni}variate medians; similar to the version \eqn{V2} of the reference;
	robust but not affine invariant.}
      \item{"random"}{based on a random selection, i.e., reproducible only
	via \code{\link{set.seed}()}.}
      \item{"manual"}{based on manual selection; in this case, a vector
	\code{man.sel} containing the indices of the selected
	observations must be specified.}
      \item{"V2"}{based on the Euclidean norm from the \bold{uni}variate
    medians; this is the version \eqn{V2} of the reference; robust but not
    affine invariant.}
    }
    \code{"Mahalanobis"} and \code{"V2"} where proposed by Hadi
    and the other authors in the reference as versions \sQuote{V_1}
    and \sQuote{V_2}, as well as \code{"manual"},
    while \code{"random"} is provided in order to study the behaviour of
    BACON.  Option \code{"dUniMedian"} is similar to \code{"V2"} and is due
    to U. Oetliker.
  }
  \item{man.sel}{only when \code{init.sel == "manual"}, the indices of
    observations determining the initial basic subset (and \code{m <-
      length(man.sel)}).}
  \item{maxsteps}{maximal number of iteration steps.}
  \item{allowSingular}{logical indicating a solution should be sought
    also when no matrix of rank \eqn{p} is found.}
  \item{verbose}{logical indicating if messages are printed which trace
    progress of the algorithm.}
}
\details{
  Remarks on the tuning parameter \code{alpha}: Let \eqn{\chi^2_p}
  be a chi-square distributed random variable with \eqn{p} degrees
  of freedom (\eqn{p} is the number of variables; \eqn{n} is the
  number of observations).  Denote the \eqn{(1-\alpha)} quantile by
  \eqn{\chi^2_p(\alpha)}, e.g., \eqn{\chi^2_p(0.05)} is the 0.95 quantile.
  Following Billor et al. (2000), the cutoff value for the
  Mahalanobis distances is defined as \eqn{\chi_p(\alpha/n)} (the square
  root of \eqn{chi^2_p}) times a correction factor \eqn{c(n,p)},
  \eqn{n} and \eqn{p},
  and they use \eqn{\alpha=0.05}.
}
\value{
  a \code{\link{list}} with components
  \item{subset}{logical vector of length \code{n} where the \code{i}-th
    entry is true iff the i-th observation is part of the final selection.}
  \item{dis}{numeric vector of length \code{n} with the (Mahalanobis)
    distances.}
  \item{cov}{\eqn{p \times p}{p x p} matrix, the corresponding robust
    estimate of covariance.}
}
\references{
  Billor, N., Hadi, A. S., and Velleman , P. F. (2000).
  BACON: Blocked Adaptive Computationally-Efficient Outlier Nominators;
  \emph{Computational Statistics and Data Analysis} \bold{34}, 279--298.
  \doi{10.1016/S0167-9473(99)00101-2}
}
\author{
  Ueli Oetliker, Swiss Federal Statistical Office, for S-plus 5.1.
  Port to \R, testing etc, by Martin Maechler;
  Init selection \code{"V2"} and correction of default \code{alpha} from 0.95 to 0.05,
  by Tobias Schoch, FHNW Olten, Switzerland.
}
\seealso{
  \code{\link[robustbase]{covMcd}} for a high-breakdown (but more computer
  intensive) method;
  \code{\link{BACON}} for a \dQuote{generalization}, notably to
  \emph{regression}.
}
\examples{
 require(robustbase) # for example data and covMcd():
 ## simple 2D example :
 plot(starsCYG, main = "starsCYG  data  (n=47)")
 B.st <- mvBACON(starsCYG)
 points(starsCYG[ ! B.st$subset,], pch = 4, col = 2, cex = 1.5)
 stopifnot(identical(which(!B.st$subset), c(7L,11L,20L,30L,34L)))
 ## finds the 4 clear outliers (and 1 "borderline");
 ## it does not find obs. 14 which is an outlier according to covMcd(.)

 iniS <- setNames(, eval(formals(mvBACON)$init.sel)) # all initialization methods, incl "random"
 set.seed(123)
 Bs.st <- lapply(iniS[iniS != "manual"], function(s)
                 mvBACON(as.matrix(starsCYG), init.sel = s, verbose=FALSE))
 ii <- - match("steps", names(Bs.st[[1]]))
 Bs.s1 <- lapply(Bs.st, `[`, ii)
 stopifnot(exprs = {
    length(Bs.s1) >= 4
    length(unique(Bs.s1)) == 1 # all 4 methods give the same
 })

 ## Example where "dUniMedian" and "V2" differ :
 data(pulpfiber, package="robustbase")
 dU.plp <- mvBACON(as.matrix(pulpfiber), init.sel = "dUniMedian")
 V2.plp <- mvBACON(as.matrix(pulpfiber), init.sel = "V2")
 (oU <- which(! dU.plp$subset))
 (o2 <- which(! V2.plp$subset))
 stopifnot(setdiff(o2, oU) \%in\% c(57L,58L,59L,62L))
 ## and 57, 58, 59, and 62 *are* outliers according to covMcd(.)

 ## 'coleman' from pkg 'robustbase'
 coleman.x <- data.matrix(coleman[, 1:6])
 Cc <- covMcd (coleman.x) # truly robust
 summary(Cc) # -> 6 outliers (1,3,10,12,17,18)
 Cb1 <- mvBACON(coleman.x) ##-> subset is all TRUE hmm??
 Cb2 <- mvBACON(coleman.x, init.sel = "dUniMedian")
 stopifnot(all.equal(Cb1, Cb2))
 ## try 20 different random starts:
 Cb.r <- lapply(1:20, function(i) { set.seed(i)
                     mvBACON(coleman.x, init.sel="random", verbose=FALSE) })
 nm <- names(Cb.r[[1]]); nm <- nm[nm != "steps"]
 all(eqC <- sapply(Cb.r[-1], function(CC) all.equal(CC[nm], Cb.r[[1]][nm]))) # TRUE
 ## --> BACON always  breaks down, i.e., does not see the outliers here
 \dontshow{stopifnot(Cb1$subset, Cb.r[[1]]$subset, eqC)}
 ## breaks down even when manually starting with all the non-outliers:
 Cb.man <- mvBACON(coleman.x, init.sel = "manual",
                   man.sel = setdiff(1:20, c(1,3,10,12,17,18)))
 which( ! Cb.man$subset) # the outliers according to mvBACON : _none_
}
\keyword{multivariate}
\keyword{robust}
