\name{lmrob.split}
\alias{lmrob.split}
\title{
  Split Continuous and Categorical Predictors
}
\description{
  Splits the design matrix into categorical and continuous
  predictors. Categorical variables are variables that are factors or
  ordered factors.
}
\usage{
lmrob.split(mf, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mf}{
    model frame (as returned by \code{model.frame})
  }
  \item{x}{
    design matrix (optional)
  }
}
\details{
  Interactions of continuous and categorical variables are
  considered as categorical variables. The M-S algorithm requires that
  continuous variables that are part of an interaction with categorical
  variables are considered categorical variables as well.
}
\value{
  A list that includes the following components:
  \item{x1 }{design matrix containing only categorical variables}
  \item{x1.idx }{logical vectors of the variables considered
    categorical in the original design matrix}
  \item{x2 }{design matrix containing the continuous variables}
}
\references{
  Maronna, R. A., and Yohai, V. J. (2000).
  Robust regression with both continuous and categorical predictors.
  \emph{Journal of Statistical Planning and Inference} \bold{89}, 197--214.
}
\author{
  Manuel Koller
}
\seealso{
  \code{\link{lmrob.M.S}}
}
\examples{
data(education)
education <- within(education, Region <- factor(Region))
obj <- lm(Y ~ Region + X1 + X2 + X3, education)
splt <- lmrob.split(obj$model)
str(splt)
}
