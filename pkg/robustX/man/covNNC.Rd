\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\name{covNNC}
\alias{covNNC}
\alias{cov.nnve}
\title{Robust Covariance Estimation via Nearest Neighbor Cleaning}
\description{
  \code{covNNC()} estimates robust covariance/dispersion matrices by the
  nearest neighbor variance estimation (NNVE) or (rather)
  \dQuote{Nearest Neighbor Cleaning} (NNC) method of Wang and Raftery
  (2002, \emph{JASA}).
}
\usage{
covNNC(X, k = min(12, n - 1), pnoise = 0.05, emconv = 0.001,
       bound = 1.5, extension = TRUE, devsm = 0.01)
}
\arguments{
  \item{X}{matrix in which each row represents an observation or point
    and each column represents a variable.}
  \item{k}{desired number of nearest neighbors (default is 12)}
  \item{pnoise}{percent of added noise}
  \item{emconv}{convergence tolerance for EM}
  \item{bound}{value used to identify surges in variance caused by
               outliers wrongly included as signal points (\code{bound = 1.5}
               means a 50 percent increase)}
  \item{extension}{whether or not to continue after reaching the last
                   chi-square distance.  The default is to continue,
                  which is indicated by setting \code{extension = TRUE}.}
  \item{devsm}{when \code{extension = TRUE}, the algorithm stops if the
      relative difference in variance is less than \code{devsm}.
      (default is 0.01)}
}

\value{
 A list with components
 \item{cov}{covariance matrix}
 \item{mu}{mean vector}
 \item{postprob}{posterior probability}
 \item{classification}{classification (0=noise otherwise 1) obtained
   by rounding \code{postprob}}
 \item{innc}{list of initial nearest neighbor cleaning results
   (components are the covariance, mean, posterior probability
   and classification)}
}
\references{
  Wang, N. and Raftery, A. (2002)
  Nearest neighbor variance estimation (NNVE):
  Robust covariance estimation via nearest neighbor cleaning (with discussion).
  \emph{Journal of the American Statistical Association} \bold{97}, 994--1019.

  see also University of Washington Statistics Technical Report 368 (2000)
  \url{http://www.stat.washington.edu/www/research/reports}
}

\author{Naisyin Wang \email{nwang@stat.tamu.edu} and
  Adrian Raftery \email{raftery@stat.washington.edu}
  with contributions from
  Chris Fraley \email{fraley@stat.washington.edu}.

  \code{covNNC()}, then named \code{cov.nnve()}, used to be (the only
  function) in CRAN package \pkg{covRobust} (2003), which was archived
  in 2012.

  Martin Maechler allowed \code{ncol(X) == 1},
  sped up the original code, by reducing the amount of scaling;
  further, the accuracy was increased (using internal \code{q.dDk()}).
  The original version is available, unexported as
  \code{robustX:::covNNC1}.% --> ../R/covNNC-orig.R
}
\seealso{
  \code{\link[MASS]{cov.mcd}} from package \CRANpkg{MASS};
  \code{\link[robustbase]{covMcd}}, and \code{\link[robustbase]{covOGK}}
  from package \CRANpkg{robustbase}.

  The whole package \CRANpkg{rrcov}.
}
\note{
  Terms of use: GPL version 2 or newer.

  MM: Even though \code{covNNC()} is backed by a serious scientific
  publication, I cannot recommend its use at all.
}

\examples{
data(iris)
covNNC(iris[-5])

data(hbk, package="robustbase")
hbk.x <- data.matrix(hbk[, 1:3])
covNNC(hbk.x)
}
\keyword{multivariate}
\keyword{robust}
