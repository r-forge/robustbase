\name{mvBACON}
\Rdversion{1.1}
\alias{mvBACON}
\title{BACON: Blocked Adaptive Computationally-Efficient Outlier Nominators}
\description{
  This function performs an outlier identification
  algorithm to the data in the x array [n x p] and y vector [n]
  following the lines described by Hadi et al. for their
  BACON outlier procedure.
}
\usage{
mvBACON(x, collect = 4, m = min(collect * p, n * 0.5), alpha = 0.95,
        init.sel = c("Mahalanobis", "dUniMedian", "random", "manual"),
        man.sel, maxsteps = 100, allowSingular = FALSE, verbose = TRUE)
}
\arguments{
  \item{x}{numeric matrix (of dimension \eqn{[n x p]}),
    not supposed to contain missing values.}
  \item{collect}{a multiplication factor, when \code{init.sel} is not
    \code{"manual"}, to define \eqn{m}, the size of the initial basic
    subset, as \code{m <- min(p * collect, n/2)}.}
  \item{m}{integer in \code{1:n} specifying the size of the initial basic
    subset; used only when \code{init.sel} is not \code{"manual"}.}
  \item{alpha}{significance level for the \eqn{\chi^2}{chisq} cutoff,
    used to define the next iterations basic subset.}
  \item{init.sel}{character string, specifying the initial selection
    mode; implemented modes are:
    \describe{
      \item{"Mahalanobis"}{based on Mahalanobis distances (default)}
      \item{"dUniMedian"}{based on the distances from the
	\bold{uni}variate medians}
      \item{"random"}{based on a random selection}
      \item{"manual"}{based on manual selection; in this case, a vector
	\code{man.sel} containing the indices of the selected
	observations must be specified.}
    }
    \code{"Mahalanobis"}, \code{"dUniMedian"} where proposed by Hadi
    and the other authors in the reference as versions \sQuote{V_1}
    and \sQuote{V_2}, as well as \code{"manual"},
    while \code{"random"} is provided in order to study the behaviour of
    BACON.
  }
  \item{man.sel}{only when \code{init.sel == "manual"}, the indices of
    observations determining the initial basic subset (and \code{m <-
      length(man.sel)}).}
  \item{maxsteps}{maximal number of iteration steps.}
  \item{allowSingular}{logical indicating a solution should be sought
    also when no matrix of rank \eqn{p} is found.}
  \item{verbose}{logical indicating if messages are printed which trace
    progress of the algorithm.}
}
% \details{
%
% }
\value{
  a list with components
  \item{subset}{logical vector of length \code{n} where the \code{i}-th
    entry is true iff the i-th observation is part of the final selection.}
  \item{dis}{numeric vector of length \code{n} with the (Mahalanobis)
    distances.}
  \item{cov}{\eqn{p \times p}{p x p} matrix, the corresponding robust
    estimate of covariance.}
}
\references{
  Billor, N., Hadi, A. S., and Velleman , P. F. (2000).
  BACON: Blocked Adaptive Computationally-Efficient Outlier Nominators;
  \emph{Computational Statistics and Data Analysis} \bold{34}, 279--298.
}
\author{
  Ueli Oetliker, Swiss Federal Statistical Office, for S-plus 5.1.
  Port to \R, testing etc, by Martin Maechler
}
% \note{
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[robustbase]{covMCD}} for a high-breakdown (but more computer
  intensive) method;
  \code{\link{BACON}} for a \dQuote{generalization}, notably to
  \emph{regression}.
}
\examples{
## simple 2D example :
 plot(starsCYG, main = "starsCYG  data  (n=47)")
 B.st <- mvBACON(starsCYG)
 points(starsCYG[ ! B.st$subset,], pch = 4, col = 2, cex = 1.5)
 ## finds the clear outliers (and 3 "borderline")

 ## 'coleman' from pkg 'robustbase'
 coleman.x <- data.matrix(coleman[, 1:6])
 Cc <- covMcd (coleman.x) # truely robust
 Cb1 <- mvBACON(coleman.x) ##-> subset is all TRUE hmm??
 Cb2 <- mvBACON(coleman.x, init.sel = "dUniMedian")
 ## --> BACON "breaks down" here
}
\keyword{multivariate}
\keyword{robust}
