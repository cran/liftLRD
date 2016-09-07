\name{artificial.levels}
\alias{artificial.levels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{artificial.levels}
\description{
     This function splits the coefficients into levels according to either (i) increasing quantiles of the removed interval lengths or (ii) dyadic splitting relative 
	to a fixed lowest scale    
}
\usage{
artificial.levels(y, rem, time, tail = TRUE, type = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a vector of the removed interval lengths (in the order of \env{removelist}). }
  \item{rem}{vector of indices of the removed points (from the output of the forward transform). }
  \item{time}{vector of observed times for the decomposed signal. }
  \item{tail}{a boolean variable indicating whether coarse levels with a small number of detail coefficients (less than 10 coefficients) should be combined.} 
  \item{type}{an integer indicating which type of artificial levels to compute.  If \code{type==1}, the "usual" artificial levels using quantiles of the scales; otherwise the levels
		are computed using quantiles being a multiplicative factor of a fixed number.  See details for more information. 
 }
}
\details{
     The function computes the so-called artificial levels of a set of removed integrals and corresponding detail coefficients, to mimic the dyadic level splitting in a classical 
     wavelet framework.  Details on the "usual" quantile-based splitting can be found in \code{\link[adlift]{artlev}}.  If \code{type==2} or \code{type==3}, the artificial levels are 
	defined by intervals of the form [a0 2^j,a0 2^(j-1) ) as described in Jansen et al. (2009), with a0 = 0.5 for \code{type==2} and set to the minimum sampling interval for 
	\code{type==3}.  The amalgamation of coarser artificial levels prevents variable energies at coarser scales affecting the predicted relationship between the wavelet scales 
	and their corresponding energies.
}
\value{
    \item{p}{a list of the grouped indices of \env{removelist} (in decreasing group size) indicating thresholding groups.}
  }

\author{Matt Nunes, Marina Knight}

\references{
Jansen, M, Nason, G. P. and Silverman, B. W. (2009) Multiscale methods for data on graphs and irregular
multidimensional situations. \emph{J. Roy. Stat. Soc. B} \bold{71}, Part 1, 97--125.
}


\seealso{\code{\link{liftHurst}} }

\examples{
#create test signal data
#
library(adlift)
x<-runif(100)
y<-make.signal2("blocks",x=x)
#
#perform forward transform...
#
out<-fwtnp(x,y,LocalPred=AdaptNeigh,neighbours=2)
#
al<-artificial.levels(out$lengthsremove,out$removelist, x, type = 1)
#
#
# the indices of removelist split into levels:
al
#
}
\keyword{arith}% at least one, from doc/KEYWORDS
%\keyword{lists}%   ONLY ONE   keyword per line









