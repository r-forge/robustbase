\name{weights.lmrob}
\title{Extract Robustness and Model Weights}
\alias{weights.lmrob}
\alias{weights.glmrob}
\description{
  \code{weights()} extracts robustness weights or fitting
  (or prior) weights from a \code{lmrob} or \code{glmrob} object.
}
\usage{
\method{weights}{lmrob}(object, type = c("prior", "robustness", "working"), ...)

\method{weights}{glmrob}(object, type = c("prior", "robustness", "working"), ...)
}
\arguments{
  \item{object}{
    an object of class \code{"lmrob"} or \code{"glmrob"}, typically the
    result of a call to \code{\link{lmrob}}, or \code{\link{glmrob}},
    respectively.}
  \item{type}{the type of weights to be returned.  Either
    \code{"prior"} (default), \code{"robustness"}, or \code{"working"};
    working weights for \code{\link{lmrob}} objects are the product of the
    prior (often all \eqn{ = 1}) and the robustness weights.
    For \code{\link{glmrob}} objects, they currently are the same as the
    \code{"robustness"} ones %% not using prior weights as long as not properly implemented there
  }
  \item{\dots}{not used currently.}
}
\details{
  The \dQuote{prior weights} correspond to the weights specified using
  the \dQuote{weights} argument when calling \code{lmrob}. The
  \dQuote{robustness weights} are the weights assigned by the
  M-estimator of regression, \eqn{\psi(r_i/S) / (r_i/S)}. The robust
  coefficient estimate then numericarlly corresponds to a weighted least
  squares fit using the product of both types of weights as weights;
  these are also called \dQuote{working weights}, and are available as
  \code{weights(<lmrob>, type = "working")}.
}
\value{
  Weights extracted from the object \code{object}.
}
\author{Manuel Koller and Martin Maechler.}

\seealso{
  \code{\link{lmrob}}, \code{\link{glmrob}} and \code{\link{weights}}
}
