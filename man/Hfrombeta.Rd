\name{Hfrombeta}
\alias{Hfrombeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Hurst exponent from wavelet scale - energy regression slope 
}
\description{
Uses the slope of the relationship between wavelet scale and wavelet energy to compute an estimate of the Hurst exponent
}
\usage{
Hfrombeta(beta, model = c("FBM","FGN","ID"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{ The estimated slope of the relationship between wavelet scale and energy.}
  \item{model}{ The assumed long-range dependence model for the time series under analysis.}
}
\details{
	There is a theoretical linear relationship growth in the (log) wavelet energy for increasing wavelet scale.  This corresponds to the decay in the autocorrelation of a 
	(long range dependent) time series being analysed, and therefore the Hurst exponent, H.  The specific relation to H is dependent to the assumed model; in particular for a 
	Fractional Brownian motion, the relationship between H and the slope is H = abs(beta - 1)/2, whereas for Fractional Gaussian noise or dth order Fractional differenced series,
	the relationship is H = (beta+1)/2.
}
\value{
 \item{H}{The Hurst exponent, computed for a specific beta and underlying model.}
}
\references{
Knight, M. I, Nason, G. P. and Nunes, M. A. (2016) A wavelet lifting approach to long-memory estimation. \emph{Stat. Comput. (to appear)} DOI 10.1007/s11222-016-9698-2.\cr\cr
Beran, J. et al. (2013)  Long-Memory Processes.  Springer. 
}
\author{
Matt Nunes 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{liftHurst}}
}
\examples{

Hfrombeta(0.8,model="FGN")

}
\keyword{manip}
