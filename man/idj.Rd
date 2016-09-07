\name{idj}
\alias{idj}
\alias{meanj}
\alias{medj}
\alias{mean2}
\alias{mad2}
\alias{meanmo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to perform summary calculations of wavelet scales and energies.
}
\description{
To estimate the slope of the relationship between wavelet scale and wavelet energy, choices have to be made as to how these quantities are computed.  Examples of these choices are 
the functions listed here.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
idj(x, j) 
meanj(x, j) 
medj(x, j) 
mean2(x) 
mad2(x) 
meanmo(x) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of values, representing the integrals or detail coefficients in a particular artificial level}
  \item{j}{the index of the artificial level.
%%     ~~Describe \code{j} here~~
}
}
\value{
A numeric value corresponding to the average squared detail coefficient, squared mean absolute deviation, median scale etc.  
}
\references{
Knight, M. I, Nason, G. P. and Nunes, M. A. (2016) A wavelet lifting approach to long-memory estimation. \emph{Stat. Comput. (to appear)} DOI 10.1007/s11222-016-9698-2.
}

\author{
Matt Nunes
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{liftHurst}}
}
\examples{

x<-rnorm(50,30,2)

# calculate the average squared value of x (i.e. energy)

mean2(x)

}
\keyword{manip}
