\name{psiFunc}
\alias{psiFunc}
\alias{huberPsi}
\alias{hampelPsi}
\title{Constructor for Objects "Psi Function" Class}
%% Rd Problem:  'psi_func' cannot be used in title!
\description{
  \code{psiFunc(..)} is a convenience interface to
  \code{new("psi_func",..)}, i.e. for constructing objects of class
  \code{"psi_func"}.
}
\usage{
psiFunc(rho, psi, wgt, Dpsi, Erho = NULL, Epsi2 = NULL, EDpsi = NULL, name, ...)

huberPsi
hampelPsi
}
\arguments{
  \item{rho, psi, wgt, Dpsi}{each a \code{\link{function}} of \code{x} and
    tuning parameters typically.}
  \item{Erho, Epsi2, EDpsi}{see \code{\link{psi_func-class}}, and note
    that these may change in the future.}
  \item{name}{Name of \eqn{\psi}{psi}-function used for printing.}
  \item{\dots}{potential further arguments for specifying tuning
    parameter names and defaults.} %% FIXME; give more details
}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
% \value{
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
% }
%%\references{ ~put references to the literature/web site here ~ }
\author{Martin Maechler}
\seealso{\code{\link{psi_func-class}} for the class description.}
\examples{
## classical {trivial, not interesting}:
F1 <- function(x, .) rep.int(1, length(x))
FF <- function(.) rep.int(1, length(.))
cPsi <- psiFunc(rho = function(x,.) x^2 / 2, psi = function(x, .) x,
                wgt = F1, Dpsi = F1,
                Erho = function(.) rep.int(1/2, length(.)),
                Epsi2 = FF, EDpsi = FF, name = "classic", . = Inf)
show(cPsi)
}
\keyword{classes}
\keyword{robust}
