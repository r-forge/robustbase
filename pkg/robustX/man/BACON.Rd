\name{BACON}
\Rdversion{1.1}
\alias{BACON}
\alias{.lmBACON}
\title{BACON for Regression or Multivariate Covariance Estimation}
\description{
  BACON, short for \sQuote{\bold{B}locked \bold{A}daptive
    \bold{C}omputationally-Efficient \bold{O}utlier \bold{N}ominators},
  is a somewhat robust algorithm (set), with an implementation for
  regression or multivariate covariance estimation.

  \code{BACON()} applies the multivariate (covariance estimation)
  algorithm, using \code{\link{mvBACON}(x)} in any case, and when
  \code{y} is not \code{NULL} adds a regression iteration phase,
  using the auxiliary \code{.lmBACON()} function.
}
\usage{
BACON(x, y = NULL, intercept = TRUE,
      m = min(collect * p, n * 0.5),
      init.sel = c("Mahalanobis", "dUniMedian", "random", "manual"),
      man.sel, init.fraction = 0, collect = 4,
      alpha = 0.95, maxsteps = 100, verbose = TRUE)

## *Auxiliary* function:
.lmBACON(x, y, intercept = TRUE,
         init.dis, init.fraction = 0, collect = 4,
         alpha = 0.95, maxsteps = 100, verbose = TRUE)
}
\arguments{
  \item{x}{a multivariate matrix of dimension [n x p]
    considered as containing no missing values.}
  \item{y}{the response (n vector) in the case of regression, or
    \code{NULL} for the multivariate case.}
  \item{intercept}{logical indicating if an intercept has to be used
    for the regression.}
  \item{m}{integer in \code{1:n} specifying the size of the initial basic
    subset; used only when \code{init.sel} is not \code{"manual"}; see
    \code{\link{mvBACON}}.}
  \item{init.sel}{character string, specifying the initial selection
    mode; see \code{\link{mvBACON}}.}
  \item{man.sel}{only when \code{init.sel == "manual"}, the indices of
    observations determining the initial basic subset (and \code{m <-
      length(man.sel)}).}
  \item{init.dis}{the distances of the x matrix used for the initial
    subset determined by \code{\link{mvBACON}}.}
  \item{init.fraction}{if this parameter is > 0 then the tedious steps
    of selecting the initial subset are skipped and an initial subset of
    size n * init.fraction is chosen (with smallest dis)}
  \item{collect}{numeric factor chosen by the user to define the size of
    the initial subset (p * collect)}
  \item{alpha}{significance level.}
  \item{maxsteps}{the maximal number of iteration steps
    (to prevent infinite loops)}
  \item{verbose}{logical indicating if messages are printed which trace
    progress of the algorithm.}
}
\details{

 init.sel: the initial selection mode; implemented modes are:
	     "Mah" -> based on Mahalanobis distance (default)
	     "dis" -> based on the distances from the medians
	     "ran" -> based on a random selection
	     "man" -> based on manual selection
		 in this case the vector 'man.sel' which contains
		 the indices of the selected observations
		 must be given.
	     "Mah" and "dis" are proposed by Hadi while "ran" and "man"
	     were implemented in order to study the behaviour of BACON.


}
\value{
  basically a \code{\link{list}} with components
  \item{subset}{the observation indices (in \code{1:n}) denoting the
    subset of ``good'' observations.}
  \item{tis}{............}
}
\references{
  Billor, N., Hadi, A. S., and Velleman , P. F. (2000).
  BACON: Blocked Adaptive Computationally-Efficient Outlier Nominators;
  \emph{Computational Statistics and Data Analysis} \bold{34}, 279--298.
}
\author{
  Ueli Oetliker, Swiss Federal Statistical Office, for S-plus 5.1;
  25.05.2001; modified six times till 17.6.2001.

  Port to \R, testing etc, by Martin Maechler.
}
\seealso{
  \code{\link{mvBACON}}, the multivariate version of the BACON
  algorithm.
}

\examples{
data(starsCYG, package = "robustbase")
## Plot simple data and fitted lines
plot(starsCYG)
 lmST <-    lm(log.light ~ log.Te, data = starsCYG)
(B.ST <- with(starsCYG,  BACON(x = log.Te, y = log.light)))
(RlmST <- lmrob(log.light ~ log.Te, data = starsCYG))
abline(lmST, col = "red")
abline(RlmST, col = "blue")

}
\keyword{robust}
\keyword{regression}
