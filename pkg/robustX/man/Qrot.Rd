\name{Qrot}
\alias{Qrot}
\title{Rotation Matrix to Specific Direction}
\description{
  Construct the \eqn{p \times p}{p x p} rotation matrix that rotates the
  unit vector (1,0,....0), i.e., the \eqn{x_1}-axis,
  onto (1,1,1,...1)/\eqn{\sqrt{p}}{sqrt(p)}, or more generally to
  \eqn{u/{\left\|u\right\|}}{u / ||u||}  (\eqn{u :=}\code{unit.image}).
}
\usage{
Qrot(p, transpose = FALSE, unit.image = rep(1, p))
}
\arguments{
  \item{p}{integer; the dimension (of the vectors involved).}
  \item{transpose}{logical indicating if the \emph{transposed} matrix is
    to returned.}
  \item{unit.image}{numeric vector of length \eqn{p} onto which the unit
    vector should be rotated; defaults to \emph{\dQuote{the diagonal}}
    \eqn{\propto}{prop. to}\eqn{(1,1,1,...,1)}.
  }
}
\details{
  The \code{\link{qr}} decomposition is used for a Gram-Schmitt basis
  orthogonalization.
}
\value{
  \eqn{p \times p}{p x p} orthogonal matrix which rotates
  \eqn{(1,0,...,0)} onto a vector proportional to \code{unit.image}.
}
\author{Martin Maechler}
\seealso{\code{\link{qr}}, matrix (and vector) multiplication,
  \code{\link{\%*\%}}.
}
\examples{
Q <- Qrot(6)
zapsmall(crossprod(Q)) # 6 x 6 unity <==> Q'Q = I  <==> Q orthogonal

if(require("MASS")) {
  Qt <- Qrot(6, transpose = TRUE)
  stopifnot(all.equal(Qt, t(Q)))
  fractions(Qt ^2) # --> 1/6  1/30  etc, in an almost lower-triagonal matrix
}
}
\keyword{array}
