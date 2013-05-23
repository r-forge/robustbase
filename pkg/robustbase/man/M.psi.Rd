\name{M.psi}
\alias{.M.chi}
\alias{.M.psi}
\alias{.M.wgt}
\alias{.M.rhoInf}
\title{Psi / Chi / Wgt / Rho Functions for *M-Estimation}
\description{
  Compute Psi / Chi / Wgt / Rho functions for M-estimation,
  i.e., including MM, etc.
  %% TODO: More, notably definitions ... but they are all nicely in the
  %% vignette.. How can we link from here to there ???

  \code{.M.rhoInf(x)} computes \eqn{\rho(\infty)}{rho(Inf)}, i.e., the
  normalizing or scaling constant for the transformation
  from \eqn{\rho(\cdot)}{rho(.)} to
  \eqn{\tilde\rho(\cdot)}{rho~(.)}, where the latter, aka as
  \eqn{\chi()}{chi()} fulfills \eqn{\tilde\rho(\infty) = 1}{rho~(Inf) = 1}
  which makes only sense for \dQuote{redescending} psi functions, i.e.,
  not for \code{"huber"}.

  \code{.M.wgt(x, *)} computes \eqn{\psi(x)/x}  (fast and numerically accurately).
}
\usage{
.M.psi(x, cc, psi, deriv = 0)
.M.chi(x, cc, psi, deriv = 0)
.M.wgt(x, cc, psi)
.M.rhoInf(cc, psi)
}
\arguments{
  \item{x}{numeric (\dQuote{abscissa} values) vector.}
  \item{cc}{numeric tuning constant, for some \code{psi} of length
    \eqn{> 1}.}
  \item{psi}{a string specifying the psi / chi / rho / wgt function;
    either \code{"huber"}, or one of the same possible specifiers as for
    \code{psi} in \code{\link{lmrob.control}}, i.e. currently,
    \code{"bisquare"}, \code{"lqq"}, \code{"welsh"}, \code{"optimal"},
    \code{"hampel"}, or \code{"ggw"}.}
  \item{deriv}{an integer, specifying the \emph{order} of derivative to
    consider; particularly, \code{.M.psi(x, *, deriv = -1)} is the
    principal function of \eqn{\psi()}{psi()}, typically denoted
    \eqn{\rho()}{rho()} in the literature.}
  %% FIXME: mention that  deriv = 2 is *partially* implemented
}
\details{
  Theoretically, as \code{.M.wgt} (defined via \code{.M.psi}),
  \code{.M.chi()} would not be needed explicitly as it can be computed
  from \code{.M.psi()} and \code{.M.rhoInf()}, namely, by
  \preformatted{.M.chi(x, *, deriv = d)  ==  .M.psi(x, *, deriv = d-1) / .M.rhoInf(*)}
  for \eqn{d = 0, 1, 2}  (and \sQuote{*} containing \code{par, psi}, and
  equality is in the sense of \code{\link{all.equal}(x,y, tol)} with a
  small \code{tol}.
}
\value{
  a numeric vector of the same length as \code{x}, with corresponding
  function (or derivative) values.
}
\references{
  See the vignette about %% \link{ .. } ????
  \dQuote{\eqn{\psi}{psi}-Functions Available in Robustbase}.
  %% ../inst/doc/psi_functions.Rnw
}
\author{
  Manuel Koller, notably for the original C implementation;
  tweaks and speedup via \code{\link{.Call}} by Martin Maechler.
}
\seealso{
 \code{\link{psiFunc}} and the \code{\linkS4class{psi_func}} class, both
 of which provide considerably more on the \R side, but are less
 optimized for speed.
}
\examples{
x <- seq(-5,7, by=1/8)
matplot(x, cbind(.M.psi(x, 4, "biweight"),
                 .M.chi(x, 4, "biweight"),
                 .M.wgt(x, 4, "biweight")), type = "l")
abline(h=0, v=0, lty=2, col=adjustcolor("gray", 0.6))

hampelPsi
(ccHa <- hampelPsi @ xtras $ tuningP $ k)
psHa <- hampelPsi@psi(x)
% FIXME: interesting as long as hampelPsi does not use .M.psi(... "hampel") !
## using .M.psi():
Mp.Ha <- .M.psi(x, cc = ccHa, psi = "hampel")
stopifnot(all.equal(Mp.Ha, psHa, tol = 1e-15))
}
\keyword{robust}
