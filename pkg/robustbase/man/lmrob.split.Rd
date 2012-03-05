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
lmrob.split(mf, x, type="f")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mf}{
    model frame (as returned by \code{model.frame})
  }
  \item{x}{
    design matrix (optional)
  }
  \item{type}{
    split type (see details)
  }
}
\details{
  Which split type is used can be controlled with the setting
  \code{split.type} in \code{\link{lmrob.control}}.

  There are three split types. The only differences between the types
  are how interactions between categorical and continuous variables are
  handled. The extra types of splitting can be used to avoid
  \emph{Too many singular resamples} errors.

  Type \code{f}, the default, assigns only the intercept, categorical and
  interactions of categorical variables to \code{x1}. Interactions of
  categorical and continuous variables are assigned to \code{x2}.

  Type \code{fi} assigns also interactions between categorical and
  continuous variables to \code{x1}.

  Type \code{fii} assigns not only interactions between categorical and
  continuous variables to \code{x1}, but also the continuous variables
  themselves. 
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
obj <- lm(Y ~ Region:X1:X2 + X1*X2, education)
lmrob.split(obj$model, type="f")$x1.idx
lmrob.split(obj$model, type="fi")$x1.idx
lmrob.split(obj$model, type="fii")$x1.idx
}
