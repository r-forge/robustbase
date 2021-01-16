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
  %%The DESCRIPTION file:
  \packageDESCRIPTION{robustX}

  %% Index:
  \packageIndices{robustX}

%% \preformatted{
%% Qrot             Rotation Matrix to Specific Direction
%% rbwheel          Multivariate Barrow Wheel Distribution
%%                  Random Vectors
%% L1median	 Compute the Multivariate L1-Median
%% BACON		 BACON for Regression or Multivariate Covariance Estimation
%% mvBACON		 BACON: Blocked Adaptive
%% 			Computationally-Efficient Outlier Nominators
%% }
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
  \url{https://cran.r-project.org/view=Robust}
}
\examples{
pairs( rbwheel(100, 4) )
}
\keyword{package}
