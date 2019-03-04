\name{mwatson.tune}
\alias{mwatson.tune}
\title{A function to calulate the bandwidth for the Watson kernel density estimator}
\description{A function to calulate the optimal bandwidth (tuning parameter) for the Watson kernel density estimator. Uses cross validation to calculate the bandwidth. ALso works for the multivariate Watson kernel density estimmator.}
\usage{mwatson.tune(x)
}
\arguments{
  \item{x}{A matrix or dataframe of values on S^(p-1). Should not contain NA's or NaN's}
  \item{low}{The lower value of the bandwdith to search.}
  \item{high}{The upper value of the bandwdith to search.}
}
\details{}
\value{
  \item{Optimal h}{The optimal bandwidth obtained through cross-validation}
  \item{cv}{The value of the maximised pseudo-likelihood}
}
\references{Michail Tsagris, Giorgos Athineou, Anamul Sajib, Eli Amson, Micah J. Waldstein (2018); Package "Directional", version 3.4}

\author{\author{Nishan Mudalige}}
\note{}
\seealso{}
\examples{
## install_github("baddstats/spherstat")
## library(spherstat)
## 
## x = rWatson(n=100, mode=c(0,0), kappa=50, win=sphwin())
## x = x[[1]]
## 
## mwatson.tune(x)
}
