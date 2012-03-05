\name{lmrob.M.S}
\alias{lmrob.M.S}
\title{
  M-S regression estimators
}
\description{
  Computes and M-S-estimator for linear regression using the
  \dQuote{M-S} algorithm.
}
\usage{
lmrob.M.S(x, y, control, mf, split)
}
\arguments{
  \item{x}{ numeric matrix for the predictors }
  \item{y}{ numeric vector for the response }
  \item{control}{ list as returned by \code{\link{lmrob.control}} }
  \item{mf}{ model.frame as returned by \code{\link{model.frame}} }
  \item{split}{ (optional) list as returned by \code{\link{lmrob.split}} }
}
\details{
  This function is used by \code{\link{lmrob}} and not intended to be
  used on its own (because an M-S-estimator has too low efficiency
  \sQuote{on its own}).

  An M-S estimator is a combination of an S-estimator for the
  continuous variables and an L1-estimator for the categorical
  variables.

  The S-estimator is estimated using a subsampling algorithm. If the
  model includes interactions between categorical and continuous
  variables, the subsampling algorithm might fail. In this case, one can
  choose to assign the interaction to the categorical side of variables
  rather than to the continuous side. This can be done via the control
  argument \code{split.type}, see \code{\link{lmrob.split}}.

  Note that the return status \code{converged} does not refer to the
  actual convergence status. The algorithm used does not guarantee
  convergence and thus true convergence is almost never reached. This
  is, however, not a problem if the estimate is only used as initial
  estimate part of an MM or SMDM estimate.
}
\value{
  A list with compontents
  \item{coefficients}{numeric vector (length \eqn{p}) of M-S-regression
    coefficient estimates.} 
  \item{scale}{the M-S-scale residual estimate}
  \item{residuals}{numeric vector (legnth \eqn{n}) of the residuals.}
  \item{weights}{numeric vector (length \eqn{n}) of the robustness weights.}
  \item{control}{the same list as the \code{control} argument.}
  \item{converged}{Convergence status (always \code{TRUE}), needed for
    \code{\link{lmrob.fit}}.}
}
\references{
  Maronna, R. A., and Yohai, V. J. (2000).
  Robust regression with both continuous and categorical predictors.
  \emph{Journal of Statistical Planning and Inference} \bold{89}, 197--214.
}
\author{
  Manuel Koller
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{lmrob}}, for a description of the available split types,
  see \code{\link{lmrob.split}}.
}
\examples{
data(education)
education <- within(education, Region <- factor(Region))
obj <- lm(Y ~ Region + X1 + X2 + X3, education)
x <- model.matrix(obj)
y <- education$Y
control <- lmrob.control()
str(lmrob.M.S(x, y, control, obj$model))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ M-S }
\keyword{ robust }
\keyword{ regression }% __ONLY ONE__ keyword per line
