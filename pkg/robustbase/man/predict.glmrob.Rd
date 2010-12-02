\name{predict.glmrob}
\alias{predict.glmrob}
\title{Predict Method for Robust GLM ("glmrob") Fits}
\description{
  Obtains predictions and optionally estimates standard errors of those
  predictions from a fitted \emph{robust} generalized linear model (GLM)
  object.
}
\usage{
\method{predict}{glmrob}(object, newdata = NULL,
       type = c("link", "response", "terms"), se.fit = FALSE,
       dispersion = NULL, terms = NULL, na.action = na.pass, \dots)
}
\arguments{
%% the following is +- copy-pasted from  predict.glm.Rd:
  \item{object}{a fitted object of class inheriting from \code{"glmrob"}.}
  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict.  If omitted, the fitted linear predictors are used.}
  \item{type}{the type of prediction required.  The default is on the
    scale of the linear predictors; the alternative \code{"response"}
    is on the scale of the response variable.  Thus for a default
    binomial model the default predictions are of log-odds (probabilities
    on logit scale) and \code{type = "response"} gives the predicted
    probabilities.  The \code{"terms"} option returns a matrix giving the
    fitted values of each term in the model formula on the linear predictor
    scale.

    The value of this argument can be abbreviated.
  }
  \item{se.fit}{logical switch indicating if standard errors are required.}
  \item{dispersion}{the dispersion of the GLM fit to be assumed in
    computing the standard errors.  If omitted, that returned by
    \code{summary} applied to the object is used.}
  \item{terms}{with \code{type="terms"} by default all terms are returned.
    A character vector specifies which terms are to be returned}
  \item{na.action}{function determining what should be done with missing
    values in \code{newdata}.  The default is to predict \code{NA}.}
  \item{\dots}{optional further arguments, currently simply passed to
    \code{\link{predict.lmrob}()}.}
}
% \details{
% If necessary, more details than the description above ~~
% }
\value{
%% the following is +- copy-pasted from  predict.glm.Rd:
  %% also correct,here ?
  If \code{se = FALSE}, a vector or matrix of predictions.
  If \code{se = TRUE}, a list with components
  \item{fit}{Predictions}
  \item{se.fit}{Estimated standard errors}
  \item{residual.scale}{A scalar giving the square root of the
    dispersion used in computing the standard errors.}
}
\author{Andreas Ruckstuhl}
\seealso{
  \code{\link{glmrob}()} fit these robust GLM models;
  \code{\link{predict.lm}()}, the method used for a non-robust fit.
}
% \examples{
% }
\keyword{models}
\keyword{regression}
