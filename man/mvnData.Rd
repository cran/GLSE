\name{mvnData}
\alias{mvnData}
\title{
Multivariate normal distribution data
}
\description{ The function provides multivariate normal distribution data for regression purposes. It can generate predictors and a response at the same time or only predictors.
}
\usage{
mvnData(sigma, n, mean, method = c("eigen", "svd", "chol"), sd = 1, Beta = NULL)
}

\arguments{
  \item{sigma}{
A variance covariance matrix of dimension \code{v} x  \code{v}.
}
  \item{n}{
The number of observations required to be generated.
}
  \item{mean}{
A vector of length \code{v}.
}
  \item{method}{
Possible methods are eigenvalue decomposition \code{"eigen"}, which is the default, singular value decomposition \code{"svd"}, and Cholesky decomposition \code{"chol"}. The method determines the matrix root of sigma.}

  \item{sd}{
A vector of standard deviations of length \code{v}, which is only needed when the response is generated.
}
  \item{Beta}{
A vector of length \code{v}, which is only needed when the response is generated.
}

}

\value{
A data frame consisting of
 \item{X }{ The generated predictors.}
 \item{Y }{ The generated response if required.}
}
\references{
 Alan, G. Frank, B. (2009), Computation of Multivariate Normal and t Probabilities.
  Lecture Notes in Statistics, Vol. 195., Springer-Verlag, Heidelberg. ISBN 978-3-642-01688-2.
  
   Alan, G. Frank, B. Tetsuhisa, M. Xuefe, M. Friedrich, L. Fabian, S. and Torsten,
  H. (2019). mvtnorm: Multivariate Normal and t Distributions. R package version 1.0-7.
  URL http://CRAN.R-project.org/package=mvtnorm.
  
    Aldahmani, S. and Dai, H. (2015). Unbiased Estimation for Linear Regression When n< v. International Journal of Statistics and Probability, 4(3), p61.

 Csardi, G., and Nepusz, T. (2006). The igraph software package for complex network research. InterJournal, Complex Systems, 1695(5), 1-9.
 
    Dethlefsen, C., and H?jsgaard, S. (2005). A common platform for graphical models in R: The gRbase package. Journal of Statistical Software, 14(17), 1-12.
    
}


\seealso{
 \code{\link{mvtData}}
}
\examples{
library(mvtnorm)
# Generate sigma

set.seed(1234)
sigma <- matrix(0, nrow =5,ncol=5)
sigma[1:5,1:5]<-.5
diag(sigma)<-3


# Generate  vectors of Beta, mean and standard deviation

Beta1<-  round(runif(5,1.5,3.5),1)
mu <-  runif(5, 0, 0)
sd <-  runif(5, 1, 1)

# Get the multivariate normal distribution data

set.seed(123)
result <-mvnData(sigma=sigma, n=10,mean=mu,sd=1,Beta=Beta1)

result

}

\keyword{ package }
