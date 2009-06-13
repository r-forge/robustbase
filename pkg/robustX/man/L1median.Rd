\name{L1median}
\title{Compute the Multivariate L1-Median}
\Rdversion{1.1}
\alias{L1median}
\alias{optimMethods}
\alias{nlminbMethods}
\description{
  Compute the multivariate \eqn{L_1}{L1}-median \eqn{m}, i.e., the
  minimizer of
  \deqn{\sum_{i=1}^n \| x_i - m \|,}{sum(i = 1:n; || x[i] - m ||),}
  where \eqn{\|u\| = \sqrt{\sum_{j=1}^p u_j^2}}{%
      || u || = sqrt(sum(j = 1:p; u[j]^2))}.

  As a convex problem, there's always a global minimizer, computable not
  by a closed formula but rather an iterative search.  As the
  (partial) first derivatives of the objective function is
  undefined the data points, the minimization is not entirely trivial.
}
\usage{
L1median(X, m.init = apply(X, 2, median), weights = NULL,
	method = c("nlm", "HoCrJo", "VardiZhang", optimMethods, nlminbMethods),
	pscale = apply(abs(centr(X, m.init)), 2, mean, trim = 0.40),
	tol = 1e-08, maxit = 200, trace = FALSE,
	zero.tol = 1e-15, \dots)
}
\arguments{
  \item{X}{numeric \code{\link{matrix}} of dimension \eqn{n \times p}{n x p},
    say.}
  \item{m.init}{starting value for \eqn{m}; typically and by default the
    coordinatewise median.}
  \item{weights}{optional numeric vector of non-negative weights;
    currently only implemented for method \code{"VardiZhang"}.}
  \item{method}{character string specifying the computational method,
    i.e., the algorithm to be used (can be abbreviated).}
  \item{pscale}{numeric p-vector of positive numbers,
    the coordinate-wise scale (typical size of
    \eqn{\delta{m_j}}{delta(m[j])}), where \eqn{m} is the problem's solution.}
  \item{tol}{positive number specifying the (relative) convergence
    tolerance.}
  \item{maxit}{positive integer specifying the maximal number of
    iterations (before the iterations are stopped prematurely if necessary).}
  \item{trace}{an integer specifying the tracing level of the
    iterations; \code{0} does no tracing}
  \item{zero.tol}{for method \code{"VardiZhang"}, a small positive
    number specifying the tolerance for determining that the iteration
    is \sQuote{exactly} at a data point (which is singularity).}
  \item{\dots}{optional arguments to \code{\link{nlm}()} or the
    \code{control} (list) arguments of \code{\link{optim}()}, or
    \code{\link{nlminb}()}, respectively.}
}
\details{
  Currently, we have to refer to the \dQuote{References} below.
  %% say more,  __ FIXME __
}
\value{
  currently the result \emph{depends} strongly on the \code{method}
  used.

  FIXME. This will change considerably.%% __ FIXME __
}
\references{
 Hossjer and Croux, C. (1995).
 Generalizing Univariate Signed Rank Statistics for Testing
 and Estimating a Multivariate Location Parameter.
 \emph{Non-parametric Statistics} \bold{4}, 293--308.

 Vardi, Y. and Zhang, C.-H. (2000).
 The multivariate \eqn{L_1}-median and associated data depth.
 \emph{Proc. National Academy of Science} \bold{97}(4), 1423--1426.
}
\author{
  Martin Maechler.  Method \code{"HoCrJo"} is mostly based on Kristel
  Joossens' function, implementing Hossjer and Croux (1995).
}
\seealso{
  \code{\link{median}}, \code{\link[robustbase]{covMcd}}
}
\examples{
data(stackloss)
L1median(stackloss)
L1median(stackloss, method = "HoCrJo")
%% more interesting stuff --> ../tests/L1med-tst.R
}
\keyword{robust}
\keyword{multivariate}