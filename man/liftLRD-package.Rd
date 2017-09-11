\name{liftLRD-package}
\alias{liftLRD-package}
\alias{liftLRD}
\docType{package}
\title{
Wavelet lifting estimators of the Hurst exponent for regularly and irregularly sampled time series
}
\description{
Implementations of Hurst exponent estimators based on the relationship between wavelet lifting scales and wavelet energy
}
\details{
Package information:
\tabular{ll}{
Package: \tab liftLRD\cr
Type: \tab Package\cr
Version: \tab 1.0-5\cr
Date: \tab 18/08/2016\cr
License: \tab GPL-2\cr
}

This package exploits a wavelet transform for irregularly spaced data to form wavelet-like scale-based energy measures for a time series.  This is then used to estimate the Hurst exponent. The main routine is \cr\cr
\code{\link{liftHurst}}
 
}
\references{
Knight, M. I, Nason, G. P. and Nunes, M. A. (2017) A wavelet lifting approach to long-memory estimation. \emph{Stat. Comput.} \bold{27} (6), 1453--1471. DOI 10.1007/s11222-016-9698-2.\cr\cr
For related literature on the lifting methodology adopted in the technique, see\cr\cr
Nunes, M. A., Knight, M. I and Nason, G. P. (2006) Adaptive lifting for nonparametric regression. \emph{Stat. Comput.} \bold{16} (2), 143--159.\cr\cr 
Knight, M. I. and Nason, G. P. (2009) A 'nondecimated' wavelet transform. \emph{Stat. Comput.} \bold{19} (1), 1--16. \cr\cr
For more information on long-memory processes, see e.g.\cr\cr
Beran, J. et al. (2013) Long-memory processes.  Springer.
}
\author{
Marina Knight, Guy Nason, Matt Nunes 

Maintainer: Matt Nunes <m.nunes@lancaster.ac.uk>
}

\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[adlift]{fwtnp}}, \code{\link[nlt]{fwtnpperm}} ~~
%}
