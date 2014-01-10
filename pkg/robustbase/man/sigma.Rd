\name{sigma}
\alias{sigma}
\title{Extract Residual Standard Error 'Sigma'}
\description{
  Extract the residual standard error from a fitted model.
}
\usage{
  sigma(object, ...)
}
\arguments{
  \item{object}{a fitted model.}

  \item{...}{additional, optional arguments.  (None are
  used in our methods)}
}
\value{
  the residual standard error as a scalar
}
\details{
  This is a (S3) generic function (as e.g., in package \pkg{lme4}).
  Currently, we provide methods for \code{\link{lmrob}} and
  \code{\link{nlrob}}.
}

