\name{lmrob.control}
\alias{lmrob.control}
\alias{lmrob.control.sfs}
\title{ Tuning parameters for lmrob }
\description{
  Tuning parameters for \code{\link{lmrob}}, the MM-type regression
  estimator and the associated S-, M- and D-estimators. The function
  \code{lmrob.control.sfs} provides the same control list, but with the
  default as suggested by Koller and Stahel (2010)
}
\usage{
lmrob.control(seed = NULL, nResample = 500,
              tuning.chi = NULL, bb = 0.5, tuning.psi = NULL,
              max.it = 50, groups = 5, n.group = 400,
              k.fast.s = 1, best.r.s = 2, k.max = 200,
              refine.tol = 1e-7, rel.tol = 1e-7, trace.lev = 0,
              compute.rd = FALSE, method = 'MM',
              psi = c('bisquare', 'ggw', 'welsh', 'optimal', 'hampel'),
              numpoints = 10, cov = '.vcov.avar1', ...)

lmrob.control.sfs(..., max.it = 500, k.max = 2000,
                  method = 'SMDM', psi = 'ggw', cov = '.vcov.w')
}
\arguments{
  \item{seed}{an integer vector, the seed to be used for random
    re-sampling used in obtaining candidates for the initial
    S-estimator; see \code{\link{.Random.seed}}.  The current value of
    \code{.Random.seed} will be preserved if \code{seed} is set;
    otherwise (by default), \code{.Random.seed} will be modified as
    usual from calls to \code{\link{runif}()}.
  }
  \item{nResample}{number of re-sampling candidates to be
    used to find the initial S-estimator.  Currently defaults to 500
    which works well in most situations (see references).}
  \item{tuning.chi}{tuning constant for the S-estimator.
    Sensible defaults are set according to \code{psi} to yield a 50\%
    breakdown estimator. See \emph{Details}.} 
  \item{bb}{expected value under the normal model of the
    \dQuote{chi} (rather \eqn{\rho (rho)}{rho}) function with tuning
    constant equal to \code{tuning.chi}.  This is used to compute the
    S-estimator.}
  \item{tuning.psi}{tuning constant for the re-descending M-estimator.
    Depending on the value of \code{psi} this constant is set to yield
    an estimator with asymptotic efficiency of 95\% for normal errors.
    See \emph{Details}.}
  \item{max.it}{integer specifying the maximum number of IRWLS iterations.}
  \item{groups}{(for the fast-S algorithm): Number of
    random subsets to use when the data set is large.}
  \item{n.group}{(for the fast-S algorithm): Size of each of the
    \code{groups} above.  Note that this must be at least \eqn{p}.}
  \item{k.fast.s}{(for the fast-S algorithm): Number of
    local improvement steps (\dQuote{\emph{I-steps}}) for each
    re-sampling candidate.}
  \item{best.r.s}{(for the fast-S algorithm): Number of
    of best candidates to be iterated further (i.e.,
    \dQuote{\emph{\bold{r}efined}}); is denoted \eqn{t} in
    Salibian-Barrera \& Yohai(2006).}
  \item{k.max}{(for the fast-S algorithm): maximal number of
    refinement steps for the \dQuote{fully} iterated best candidates.}
  \item{refine.tol}{(for the fast-S algorithm): relative convergence
    tolerance for the fully iterated best candidates.}
  \item{rel.tol}{(for the RWLS iterations of the MM algorithm): relative
    convergence tolerance for the parameter vector.}
  \item{trace.lev}{integer indicating if the progress of the MM-algorithm
    should be traced (increasingly); default \code{trace.lev = 0} does
    no tracing.}
  %% NOTE that lmrob.S() has its "own"  'trace.lev' !
  \item{compute.rd}{logical indicating if robust distances (based on
    the MCD robust covariance estimator \code{\link{covMcd}}) are to be
    computed for the robust diagnostic plots.  This may take some
    time to finish, particularly for large data sets, and can lead to
    singularity problems when there are \code{\link{factor}} explanatory
    variables (with many levels, or levels with \dQuote{few}
    observations). Hence, is \code{FALSE} by default.}
  \item{method}{string specifying the estimator-chain. \code{MM}
    is interpreted as \code{SM}. See \emph{Details} of
    \code{\link{lmrob}} for a description of the possible values.}
  \item{psi}{string specifying the type \eqn{\psi}-function
    used. See \emph{Details} of \code{\link{lmrob}}.}
  \item{numpoints}{Number of points used in Gauss quadrature.}
  \item{cov}{Function or string with function name to be used to
    calculate covariance matrix estimate. See \emph{Details} of
    \code{\link{lmrob}}.}
  \item{...}{Further arguments are added to the control list.}
}
\details{The default values set for \code{tuning.chi} and
  \code{tuning.psi} are:
  \tabular{rll}{
    \code{psi} \tab \code{tuning.chi} \tab \code{tuning.psi} \cr
    \code{bisquare} \tab \code{1.54764} \tab \code{4.685061} \cr
    \code{welsh} \tab \code{0.5773502} \tab \code{2.11} \cr
    \code{ggw} \tab \code{5.5} \tab \code{5.95} \cr
    \code{optimal} \tab \code{0.4047} \tab \code{1.060158} \cr
    \code{hampel} \tab \code{0.2119163} \tab \code{0.9014}
  }
}\references{
  Koller, M. and Stahel, W.A. (2010, tbp)
}
\author{ Matias Salibian-Barrera, Martin Maechler and Manuel Koller}
\seealso{ \code{\link{lmrob}}, also for references and examples.
}
\examples{
## Show the default settings:
str(lmrob.control())
}
\keyword{robust}
\keyword{regression}

