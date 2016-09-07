\name{bootci}
\alias{bootci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bootstrap confidence interval calculation
}
\description{
This function uses the Hurst exponent estimates from different lifting trajectories to form a bootstrap confidence interval 
}
\usage{
bootci(x, level = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of estimates for the Hurst exponent.
}
  \item{level}{the significance level for which to compute the bootstrap confidence interval.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A vector of length 2, indicating the lower and upper confidence interval values.
}
\references{
Knight, M. I, Nason, G. P. and Nunes, M. A. (2016) A wavelet lifting approach to long-memory estimation. \emph{Stat. Comput. (to appear)} DOI 10.1007/s11222-016-9698-2.\cr\cr
}


\author{
Matt Nunes
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{liftHurst}}
}

\examples{

library(fArma)

x<-rnorm(100,0.7,0.24)	# vector representing Hurst estimates from NLT

bootci(x)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
