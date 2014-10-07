\name{reclas}
\title{Recursive Robust Median-like Location and Scale}
\alias{reclas}
\alias{plot.reclas}
\description{
 Calculate an estimate of location, asymptotically
  equivalent to the median, and an estimate of scale
  equal to the \bold{MEAN} absolute deviation.  Both
  done recursively.
}
\usage{%-- ../R/reclas.R
reclas(y, b = 0.2, mfn = function(n) 0.1 * n^(-0.25),
     nstart = 30, m0 = median(y0),
     scon=NULL, updateScale = is.null(scon))

%% FIXME -- add these "surprising arguments":
%% \S3method{plot}{reclas}(x, M = tail(x$locn, 1), ylim = NULL,
%%      s.y = tail(x$scale, 1), se = TRUE, col.se = adjustcolor("skyblue4", 3/4),
%%      ylab = "locn", main = deparse(x$call)[1], ...)
}
\arguments{
  \item{y}{numeric vector of i.i.d. data whose location and scale
    parameters are to be estimated.}
  \item{b}{numeric tuning parameter (default value equal to that used by
    Holst, 1987).}
  \item{mfn}{a \code{\link{function}} of the index of the data which
    must be positive and and tend to 0 as the index tends to infinity.
    The default function is that used by Holst, 1987.}
  \item{nstart}{number of starting values: Starting values for the
    algorithm are formed from the first \code{nstart} values of
    \code{y}.  The default value is that used in Cameron and Turner, 1993.}
  \item{m0}{value for the initial approximate median; by default, the
    \code{\link{median}} of the first \code{nstart} observations.}
  \item{scon}{value for the scale parameter \code{s}, a
    \code{\link{function}} or \code{NULL}. When \code{NULL}, as by
    default, the scale is initialized to the mean of the absolute
    differences between the first \code{nstart} y values and \code{m0}.
    If scon is a function, the initial scale is set to \code{scon(y0, m0)},
    where \code{y0} is the vector of the first \code{nstart} y values.
    Note that \code{scon} also determines the default for \code{updateScale}.}
  \item{updateScale}{a logical indicating if the scale, initialized
    from \code{scon} should be updated in each iteration.  Otherwise, the
    the scale is held constant throughout and the algorithm becomes
    equivalent to the algorithm of Holst.}
}
\value{
  An S3 \dQuote{object} of \code{\link{class}} \code{"reclas"}; simply a
  list with entries
  \item{locn}{the successive recursive estimates of location.  The
    first \code{nstart - 1} of these are \code{NA}.}
  \item{scale}{the successive recursive estimates of scale \bold{if}
    \code{updateScale} is true; otherwise the constant value used for
    the scale.}
  \item{updateScale}{the same as the function argument.}
  \item{call}{the function call, i.e., \code{\link{match.call}}.}

  There is a \code{\link{plot}} method for \code{"reclas"}, see the
  examples.
}
\references{
  Cameron, Murray A. and Turner, T. Rolf (1993).
  Recursive location and scale estimators.
  \emph{Commun. Statist. --- Theory Meth.} \bold{22}(9) 2503--2515.

  Holst, U. (1987).
  Recursive estimators of location.
  \emph{Commun. Statist. --- Theory Meth.} \bold{16} (8) 2201--2226.

%% MM: More recent litterature -- all from computer science --
%% --> ~/R/MM/STATISTICS/robust/approx-median-large
}
\author{
 \email{r.turner@auckland.ac.nz}
 \url{http://www.stat.auckland.ac.nz/~rolf}

 Extensions by Martin Maechler (\code{scon} as function;
 \code{updateScale}, \code{plot()}).
}
\examples{
set.seed(42)
y <- rt(10000, df = 1.5) # not quite Gaussian ...
z1 <- reclas(y)
z3 <- reclas(y, scon= 1  ) # correct fixed scale
z4 <- reclas(y, scon= 100) # wrong fixed scale
z2 <- reclas(y, # a more robust initial scale:
           scon = function(y0, m0) robustbase::Qn(y0 - m0),
           updateScale = TRUE) # still updated

## Visualizing -- using the plot() method for "reclas":
M <- median(y) ; yl <- c(-1,1)* 0.5
OP <- par(mfrow=c(2,2), mar=.1+c(3,3,1,1), mgp=c(1.5, .6, 0))
 plot(z1, M=M, ylim=yl)
 plot(z2, M=M, ylim=yl)
 plot(z3, M=M, ylim=yl)
 plot(z4, M=M, ylim=yl)
par(OP)
}
\keyword{univar}
\keyword{robust}
