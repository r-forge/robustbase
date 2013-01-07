\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,psi_func-method}
\title{Plot of an object of psi_func-class}
\description{
  A simple function to visualize \code{\link{psi_func-class}}
  objects.
}
\usage{
\S4method{plot}{psi_func}(x, y,
  which = c("rho", "psi", "Dpsi", "wgt", "Dwgt"),
  shortMain = FALSE,
  col = c("black", "red3", "blue3", "dark green", "light green"),
  leg.loc = "right", ...)
}
\arguments{
  \item{x}{\code{\link{psi_func-class}} object to be plotted}
  \item{y}{(optional) vector of points to plot object at}
  \item{which}{slots to be included in plot, by default
    all of the slots are included}
  \item{shortMain}{use short or long plot title}
  \item{col}{colors to be used for the different slots}
  \item{leg.loc}{legend placement, see also \code{x} argument
    of \code{\link{legend}}}
  \item{...}{passed to \code{\link{matplot}}}
}
%% \section{Methods}{
%% \describe{
%%   \item{\code{signature(x = "psi_func", y = "ANY")}}{
%%     Plot function for \code{\link{psi_func-class}} objects.
%% }
%% }}
\examples{
plot(huberPsi)
plot(huberPsi, which=c("rho", "psi", "Dwgt"), shortMain=TRUE)
}
\keyword{methods}
