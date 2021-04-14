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
      init.sel = c("Mahalanobis", "dUniMedian", "random", "manual", "V2"),
      man.sel, init.fraction = 0, collect = 4,
      alpha = 0.05, alphaLM = alpha, maxsteps = 100, verbose = TRUE)

## *Auxiliary* function:
.lmBACON(x, y, intercept = TRUE,
         init.dis, init.fraction = 0, collect = 4,
         alpha = 0.05, maxsteps = 100, verbose = TRUE)
}
\arguments{
  \item{x}{a multivariate matrix of dimension [n x p]
    considered as containing no missing values.}
  \item{y}{the response (n vector) in the case of regression, or
    \code{NULL} for the multivariate case, where just
    \code{\link{mvBACON}()} is returned.}
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
  \item{alpha}{number in \eqn{(0, 1)} determining the cutoff value for the
    Mahalanobis distances (multivariate outlier nomination in
    \code{\link{mvBACON}()}), or the discrepancies for regression, see \code{alphaLM}.}
  \item{alphaLM}{number in \eqn{(0, 1)} where a \code{1-alphaM} t-quantile
    is the cutoff for the discrepancies (for regression, \code{.lmBACON()}); see details.}
  \item{maxsteps}{the maximal number of iteration steps
    (to prevent infinite loops)}
  \item{verbose}{logical indicating if messages are printed which trace
    progress of the algorithm.}
}
\note{\dQuote{BACON} was also chosen in honor of Francis Bacon:

  \emph{Whoever knows the ways of Nature will more easily notice her deviations;
    and, on the other hand, whoever knows her deviations will more accurately
    describe her ways.}
  \cr Francis Bacon (1620), Novum Organum II 29.
}
\details{
  Notably about the initial selection mode, \code{init.sel}, see its
  description in the \code{\link{mvBACON}} arguments list.

  The choice of \code{alpha} and \code{alphaLM}:
  \itemize{
    \item Multivariate outlier nomination: see the Details section of
    \code{\link{mvBACON}}.

    \item Regression: Let \eqn{t_r(\alpha)}{t_r(ap)} denote the \eqn{1-\alpha}{1-ap} quantile of the
    Student \eqn{t}-distribution with \eqn{r} degrees of freedom,
    where \eqn{r} is the number of elements in the current subset;
    e.g., \eqn{t_r(0.05)} is the 0.95 quantile.  Following Billor et
    al. (2000), the cutoff value for the discrepancies is defined
    as \eqn{t_r(\alpha/(2r + 2))}{t_r(ap/(2r + 2))}, and
    they use \eqn{\alpha=0.05}{ap=0.05}.   Note that this is argument
    \code{alphaLM} (defualting to \code{alpha}) for \code{BACON()}.
  }
}
\value{
  \code{BACON(x,y,..)} (for regression) returns a \code{\link{list}} with
  components
  \item{subset}{the observation indices (in \code{1:n}) denoting a
    subset of \dQuote{good} supposedly outlier-free observations.}
  \item{tis}{the \eqn{t_i(y_m, X_m)}{t[i](y[m],X[m])} of eq (6) in the
    reference; the clean \dQuote{basic subset} in the algorithm is defined
    the observations \eqn{i} with the smallest \eqn{|t_i|}{|t[i]|}, and the
    \eqn{t_i}{t[i]} can be regarded as scaled predicted errors.}
  \item{mv.dis}{the (final) discrepancies or distances of
    \code{mvBACON()}.}
  \item{mv.subset}{the \dQuote{good} subset from \code{\link{mvBACON}()},
    used to start the regression iterations.}
}
\references{
  Billor, N., Hadi, A. S., and Velleman , P. F. (2000).
  BACON: Blocked Adaptive Computationally-Efficient Outlier Nominators;
  \emph{Computational Statistics and Data Analysis} \bold{34}, 279--298.
  \doi{10.1016/S0167-9473(99)00101-2}
}
\author{
  Ueli Oetliker, Swiss Federal Statistical Office, for S-plus 5.1;
  25.05.2001; modified six times till 17.6.2001.

  Port to \R, testing etc, by Martin Maechler.
  Daniel Weeks (at pitt.edu) proposed a fix to a long standing buglet in
  \code{GiveTis()} computing the \eqn{t_i}{t[i]}, which was further improved
  Maechler, for \CRANpkg{robustX} version 1.2-3 (Feb. 2019).

  Correction of \code{alpha} default, from 0.95 to 0.05, by Tobias Schoch,
  see \code{\link{mvBACON}}.
}
\seealso{
  \code{\link{mvBACON}}, the multivariate version of the BACON
  algorithm.
}

\examples{
data(starsCYG, package = "robustbase")
## Plot simple data and fitted lines
plot(starsCYG)
lmST <- lm(log.light ~ log.Te, data = starsCYG)
abline(lmST, col = "gray") # least squares line
str(B.ST <- with(starsCYG,  BACON(x = log.Te, y = log.light)))
## 'subset': A good set of of points (to determine regression):
colB <- adjustcolor(2, 1/2)
points(log.light ~ log.Te, data = starsCYG, subset = B.ST$subset,
       pch = 19, cex = 1.5, col = colB)
## A BACON-derived line:
lmB <- lm(log.light ~ log.Te, data = starsCYG, subset = B.ST$subset)
abline(lmB, col = colB, lwd = 2)

require(robustbase)
(RlmST <- lmrob(log.light ~ log.Te, data = starsCYG))
abline(RlmST, col = "blue")
}
\keyword{robust}
\keyword{regression}
