\name{w.kde}
\alias{w.kde}
\title{A function to obtain kernel density values using the Watson kernel density estimator}
\description{A function which provides obtain kernel density estimates using the Watson kernel density estimator for a given matrix.}
\usage{
w.kde(x)
}
\arguments{
  \item{x}{A matrix or dataframe of values on S^(p-1). Should not contain NA's or NaN's}
  \item{h}{The bandwidth. If left as NULL, the optimal bandwidth will be calculated using the Watson kernel density estimator.}
}
\details{}
\value{
  \item{h}{The bandwidth.}
  \item{f}{A vector of kernel density estimates.}
}
\references{Michail Tsagris, Giorgos Athineou, Anamul Sajib, Eli Amson, Micah J. Waldstein (2018); Package "Directional", version 3.4}
\author{Nishan Mudalige}
\note{}
\seealso{}
\examples{## install_github("baddstats/spherstat")
## library(spherstat)
## 
## x = rWatson(n=100, mode=c(0,0), kappa=50, win=sphwin())
## x = x[[1]]
## 
## w.kde(x)
}
