\name{robustX-package}
\alias{robustX-package}
\alias{robustX}
\docType{package}
\title{eXperimental eXtraneous ... Functionality for Robust Statistics}
% .. : eXotic eXperts
\description{
  The package \pkg{robustX} aims to be a collection of \R functionality
  for robust statistics of methods and ideas that are considered
  as proposals, experimental, for experiences or just too much
  specialized to be part of the \dQuote{Robust Basics} package
  \pkg{robustbase}.
}
\details{
 \tabular{ll}{
  Package: \tab robustX\cr
  Date: \tab 2008-12-17, was the first, R-forge only  \cr
  Version: \tab 1.2-0\cr
  Depends: \tab robustbase\cr
  License: \tab GPL (>= 2)\cr
 }

Index:
\preformatted{
Qrot             Rotation Matrix to Specific Direction
rbwheel          Multivariate Barrow Wheel Distribution
                 Random Vectors
L1median	 Compute the Multivariate L1-Median
BACON		 BACON for Regression or Multivariate Covariance Estimation
mvBACON		 BACON: Blocked Adaptive
			Computationally-Efficient Outlier Nominators
}
}
\author{
Werner Stahel, Martin Maechler and potentially others

Maintainer: Martin Maechler
}
% \references{
% }
\seealso{
  Package \pkg{robustbase} which it complements and on which it depends;
  further package \pkg{robust} and the whole CRAN task view on robust
  statistics,
  \url{http://cran.CH.r-project.org/web/views/Robust.html}
}
\examples{
pairs( rbwheel(100, 4) )
}
\keyword{package}
