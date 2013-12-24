\name{nlrob-algorithms}
\alias{nlrob.algorithms}
\alias{nlrob.MM}
\alias{nlrob.tau}
\alias{nlrob.CM}
\alias{nlrob.mtl}
\title{MM-, Tau-, CM-, and MTL- Estimators for Nonlinear Robust Regression}
\description{
  \describe{
    \item{"MM":}{Compute an MM-estimator for nonlinear robust (constrained)
      regression.}
    \item{"tau":}{Compute a Tau-estimator for nonlinear robust (constrained)
      regression.}
    \item{"CM":}{Compute a \dQuote{Constrained M} (=: CM) estimator for
      nonlinear robust (constrained) regression.}
    \item{"MTL":}{Compute a \dQuote{Maximum Trimmed Likelihood} (=: MTL)
      estimator for nonlinear robust (constrained) regression.}
  }
}
\usage{
nlrob.MM(formula, data, pnames, lower, upper,
	 tol = 1e-06,
	 psi = c("bisquare", "lqq", "optimal", "hampel"),
         init = c("S", "lts"),
	 ctrl = nlrob.control("MM", psi = psi, init = init, fnscale = NULL,
		       tuning.chi.scale = .psi.conv.cc(psi, .Mchi.tuning.defaults[[psi]]),
		       tuning.psi.M     = .psi.conv.cc(psi, .Mpsi.tuning.defaults[[psi]]),
		       optim.control = list(), jOptimArgs = list(...)),
	 ...)

nlrob.tau(formula, data, pnames, lower, upper,
	  tol = 1e-06, psi = c("bisquare", "optimal"),
	  ctrl = nlrob.control("tau", psi = psi, fnscale = NULL,
			tuning.chi.scale = NULL, tuning.chi.tau = NULL,
			jOptimArgs = list(...)),
	  ...)

nlrob.CM(formula, data, pnames, lower, upper,
	 tol = 1e-06,
	 psi = c("bisquare", "lqq", "welsh", "optimal", "hampel", "ggw"),
	 ctrl = nlrob.control("CM", psi = psi, fnscale = NULL,
                        tuning.chi = NULL, jOptimArgs = list(...)),
	 ...)

nlrob.mtl(formula, data, pnames, lower, upper,
	  tol = 1e-06,
	  ctrl = nlrob.control("mtl", cutoff = 2.5, jOptimArgs = list(...)),
	  ...)
}
\arguments{
  \item{formula}{nonlinear regression \code{\link{formula}}, using both
    variable names from \code{data} and parameter names from \code{pnames}.}
  \item{data}{data to be used, a \code{\link{data.frame}}}
  \item{pnames}{a \code{\link{character}} vector of parameter names
    (used in \code{formula} above).}
  \item{lower, upper}{bounds aka \dQuote{box constraints} for all the
    parameters.}

  \item{tol}{numerical convergence tolerance.}
  \item{psi, init}{see \code{\link{nlrob.control}}.}
  \item{ctrl}{a \code{\link{list}}, typically the result of a call to
    \code{\link{nlrob.control}}.}

  \item{tuning.psi.M}{..}% FIXME
  \item{optim.control}{..}% FIXME
  \item{jOptimArgs}{a \code{\link{list}} of optional arguments for
    optimization, e.g., \code{trace = TRUE}, passed to to the optimizer,
    which currently must be \code{\link{JDEoptim}(.)}.}
  \item{...}{alternative way to pass the \code{jOptimArgs} above.}
}
\value{
  an \R object of \code{\link{class}} \code{"nlrob.<meth>"}, basically a
  list with components
  %% FIXME
}
\details{
  Copyright 2013, Eduardo L. T. Conceicao.  Available under
  the GPL (>= 2)
}
\author{
  Eduardo L. T. Conceicao; tweaks by Martin Maechler.
}
\source{
  For \code{"MTL"}:
  Maronna, Ricardo A., Martin, R. Douglas, and Yohai, Victor J. (2006).
  \emph{Robust Statistics: Theory and Methods} Wiley, Chichester, p. 133.
}
\references{
  \describe{
    \item{"MM":}{
      Yohai, V.J. (1987)
      High breakdown-point and high efficiency robust estimates for
      regression.
      \emph{The Annals of Statistics} \bold{15}, 642--656.
    }
    \item{"tau":}{
      Yohai, V.J., and Zamar, R.H. (1988).
      High breakdown-point estimates of regression by means of the
      minimization of an efficient scale.
      \emph{Journal of the American Statistical Association} \bold{83},
      406--413.
    }
    \item{"CM":}{
      Mendes, B.V.M., and Tyler, D.E. (1996)
      Constrained M-estimation for regression.
      In: \emph{Robust Statistics, Data Analysis and Computer Intensive
	Methods}, Lecture Notes in Statistics 109, Springer, New York, 299--320.
    }
    \item{"MTL":}{
      Hadi, Ali S., and Luceno, Alberto (1997).
      Maximum trimmed likelihood estimators: a unified approach,
      examples, and algorithms.
      Computational Statistics & Data Analysis \bold{25}, 251--272.

      Gervini, Daniel, and Yohai, Victor J. (2002).
      A class of robust and fully efficient regression estimators.
      The Annals of Statistics \bold{30}, 583--616.
    }

  }%describe
}
\examples{%% for more, --> ../tests/nlregrob-tst.R
fMM <- robustbase:::nlrob.MM( density ~ Asym/(1 + exp(( xmid -log(conc) )/scal )),
          data = DNase[DNase$Run == 1, ], pnames = c("Asym", "xmid", "scal"),
          lower = 0, upper = 3,
          ## call to nlrob.control to pass 'optim.control':
          ctrl = nlrob.control("MM", optim.control = list(trace = 1)),
                               jOptimArgs = list(trace = TRUE))

ftau <- robustbase:::nlrob.tau( density ~ Asym/(1 + exp(( xmid -log(conc) )/scal )),
           data = DNase[DNase$Run == 1, ], pnames = c("Asym", "xmid", "scal"),
           lower = 0, upper = 3,  trace = TRUE )

fCM <- robustbase:::nlrob.CM( density ~ Asym/(1 + exp(( xmid -log(conc) )/scal )),
           data = DNase[DNase$Run == 1, ], pnames = c("Asym", "xmid", "scal", "sigma"),
           lower = 0, upper = c(3, 3, 3, 0.1), trace = TRUE )

fmtl <- robustbase:::nlrob.mtl( density ~ Asym/(1 + exp(( xmid -log(conc) )/scal )),
           data = DNase[DNase$Run == 1, ], pnames = c("Asym", "xmid", "scal", "sigma"),
           lower = 0, upper = c(3, 3, 3, 0.1), trace = TRUE )

mods <- list(MM=fMM, tau=ftau, CM=fCM, MTL=fmtl)
(sts <- sapply(mods, `[[`, "status"))
stopifnot(sts == "converged")

sapply(mods, `[[`, "data") # currently 'language' %% FIXME

sapply(mods, `[[`, "coefficients")

## Compare with "traditional"  M-estimate:
fM <- nlrob(formula(fMM), data=eval(fMM$data), start = coef(fMM),
            psi = .Mwgt.psi1("bisquare"))
rbind(M=coef(fM), MM=coef(fMM)) # "similar" ..
}
