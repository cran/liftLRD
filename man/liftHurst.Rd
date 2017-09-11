\name{liftHurst}
\alias{liftHurst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Performs (non-decimated) lifting based estimation of the Hurst exponent
%%  ~~function to do ... ~~
}
\description{
	The function exploits the linear relationship in wavelet energy per scale to estimate the long range dependence parameter of a irregular time series. 
}
\usage{
liftHurst(x, grid = 1:length(x), model = "FGN", ntraj = 50, 
tradonly = FALSE, cutoffs = 0, cut.fine = TRUE, efun = meanmo, 
afun = idj, altype = 1, tail = TRUE, normalise = TRUE, 
level = 0.05, bc = TRUE, vc = TRUE, jsc = TRUE, BHonly=TRUE, 
verbose = FALSE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A time series, possibly irregularly spaced. 
}
  \item{grid}{The time samples corresponding to the time series \code{x}.}
  \item{model}{The underlying model the time series is assumed to follow.  Possible values are \code{"FBM"}, \code{"FGN"}, and \code{"ID"}.}
  \item{ntraj}{The number of lifting trajectories (bootstraps) used for the estimation of the Hurst exponent.}
  \item{tradonly}{A boolean variable, indicating whether the "traditional" minimum integral lifting trajectory should be used for the estimation (\code{tradonly=TRUE}).  In this case,
	the argument \code{ntraj} is not used.
}
  \item{cutoffs}{
	A vector indicating if the coarsest scales are to be removed when estimating the linear slope in the energy-scale relationship; for example if \code{cutoffs = 0}, 
	all levels will be used in the estimation.  For multiple cutoffs, all slopes will be calculated.
	}

  \item{cut.fine}{
	A boolean variable indicating if the finest (default) or coarsest scales are to be removed when estimating the linear slope in the energy-scale relationship; see \code{cutoffs} description for more information.
	}
  \item{efun}{A function indicating which measure of wavelet energy to use.  The default is \code{meanmo}, corresponding to the average squared detail coefficient within a particular 
	scale, scaled by (n-1).  Other measures could be \code{mad2} (squared median absolute deviation of the detail coefficients), but the average squared detail coefficient has been shown to be 
	more reliable than other measures.
}
  \item{afun}{A function indicating which measure of wavelet integral (scale) to use (corresponding to the x-axis of the energy-scale relationship).  The default is \code{idj}, 
	corresponding to the log2 of the dyadic wavelet scale, i.e an integer representing the scale.  Other measures could be \code{meanj} or \code{medj}, representing the  
	average (respectively median) integral value for each artificial scale.
}
  \item{altype}{An integer indicating which type of artificial levels to be used, see \cr 
\code{\link{artificial.levels}}} for more information.
  \item{tail}{A boolean variable indicating whether to amalgamate artificial levels with few coefficients, see \code{\link{artificial.levels}}} for more information.
  \item{normalise}{A boolean variable indicating whether to normalise the detail coefficients by their individual (induced) standard deviations (computed using the diagonal of the 
	lifting matrix). 
}
  \item{level}{
The significance level for the bootstrap confidence interval of the Hurst exponent estimate. 
}
  \item{bc}{
A boolean variable indicating whether bias-correction should take place or not, using the approximate chi^2 distribution of the wavelet energies. 
}
  \item{vc}{
A boolean variable indicating whether a weighted linear regression should be used when estimating the Hurst exponent.  If \code{TRUE}, the inverse of the variances of the approximate chi^2 distribution of the wavelet energies are used in the regression. 
}
  \item{jsc}{A boolean variable indicating whether the slope of the log-linear relationship between the artificial scales and the log of the integrals should be computed and used to reweight the estimate of the Hurst exponent. 
}
  \item{BHonly}{A boolean variable indicating whether only the Hurst estimate should be returned.  If \code{FALSE}, then the (log) energies used in the calculation is also returned.  Note that if bootstrapping is performed (\code{tradonly=FALSE}), then the energies for the last trajectory will be returned.
}
  \item{verbose}{A boolean variable indicating whether informative status messages should be printed to the screen. 
}
  \item{\dots}{Any other optional arguments to the function \code{fwtnpperm} function from the \code{nlt} package, which performs wavelet lifting of the time series \code{x} according
	to a specified lifting trajectory.}
}
\details{Wavelet lifting is performed on a time series to convert it into a set of wavelet coefficients and corresponding lifting integrals, specific to when the data were "lifted" 
during the decomposition.  The coefficients are then grouped into artificial levels, using the integrals to mimic the support of the wavelets in the classical wavelet setting, and 
therefore producing a notion of scale.  The coefficients in each artificial level are then used to compute values of the wavelet energy for a particular level.  The (slope of the)
linear relationship between the scales and their energies is then used in computing an estimate of the Hurst exponent for the series.  This procedure can be performed for multiple
(random) lifting trajectories, each producing a slightly different estimate.  
}
\value{
If \code{tradonly=TRUE}, the function returns a matrix of dimension \code{length(cutoffs)} x 2.  The first column are the slopes of the regression fits for each cutoff, whereas the 
second column are the corresponding estimates for the Hurst exponent. \cr\cr
If \code{tradonly=FALSE}, the function returns a matrix of dimension \code{length(cutoffs)} x 5.  The first column are the slopes of the regression fits for each cutoff, where the 
average is taken over the \code{ntraj} randomly generated lifting trajectories.  Similarly, the second column represents the average Hurst exponent for the cutoffs over all 
lifting paths.  The third column is the standard deviation of the \code{ntraj} Hurst estimates through performing non-decimated lifting.  The fourth and fifth columns are the lower and upper values of the bootstrap confidence interval of the Hurst exponent estimate.

If \code{BHonly=FALSE}, the routine also returns the energies and scales (on a log scale) which are used in the regression to estimate the decay properties of the spectrum (for the last lifting trajectory), as well as the weights used in the regression (if \code{vc = TRUE}).  if \code{jsc = TRUE}, the slope of the integral log-linear relationship is also returned.
}
\references{
Knight, M. I, Nason, G. P. and Nunes, M. A. (2017) A wavelet lifting approach to long-memory estimation. \emph{Stat. Comput.}  \bold{27} (6), 1453--1471. DOI 10.1007/s11222-016-9698-2.\cr\cr

For more details on the weighted linear regression and bias calculations, see e.g.\cr \cr

Veitch, D. and Abry, P. (1999) A Wavelet-Based Joint Estimator of the Parameters of Long-Range Dependence. \emph{IEEE Trans. Info. Theory} \bold{45} (3), 878--897.
}
\author{
Marina Knight, Matt Nunes
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{artificial.levels}}, \code{\link{Hfrombeta}}, \code{\link[nlt]{fwtnpperm}}
}
\examples{

library(fracdiff)

# simulate a long range dependent time series
x<-fracdiff.sim(n = 200, d = 0.3)$series

# perform lifting-based estimation of the Hurst exponent

Hestx<-liftHurst(x, tradonly=TRUE)

# An example with missingness
x1<-fracdiff.sim(n = 500, d = 0.3)$series
gap1<-30:40
gap2<-77:90
gap3<-146:166
timeindex<-setdiff(1:500,c(gap1,gap2,gap3))

Hestx1<-liftHurst(x1[timeindex],grid=timeindex, tradonly=TRUE)


}

\keyword{manip}
\keyword{regression}
\keyword{ts}
