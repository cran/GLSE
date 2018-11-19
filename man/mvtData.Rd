\name{mvtData}
\alias{mvtData}
\title{
Multivariate t-distribution data
}
\description{ The function provides multivariate t-distribution data for regression purposes. It can generate predictors and a response at the same time or only predictors.
}
\usage{
mvtData(sigma = sigma, n, mean, df = 1, delta = rep(0, nrow(sigma)),
type = c("shifted", "Kshirsagar"), sd, Beta = NULL)
}

\arguments{
  \item{sigma}{
A scale matrix of \code{v} x  \code{v} dimension.
}
  \item{n}{
The number of observations required to be generated.
}
  \item{mean}{
A vector of length \code{v}.
}
  \item{df}{
Degrees of freedom, where the default is 1.}
 \item{delta}{
A vector of noncentral parameters of length n. }

  \item{type}{
Possible types are "Kshirsagar" and "shifted". More information is provided in the R package 'mvtnorm' version 1.0-7 under the rmvt function. }


  \item{sd}{
A vector of standard deviations, which is only needed when the response is generated.
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
 \code{\link{mvnData}}
}
\examples{
library(mvtnorm)

  # Generate sigma

    sigma <- matrix(0, nrow =5,ncol=5)
    sigma[1:5,1:5]<-0
    diag(sigma)<-1


  # Generate  vectors of Beta, mean and standard deviation

    Beta1<-  round(runif(5,1.5,3.5),1)
    mu <-  runif(5, 0, 0)
    sd <-  runif(5, 1, 1)

  # Get the multivariate normal distribution data

    set.seed(123)
    result <-mvtData(sigma=sigma, n=10,mean=mu,df=1,sd=1,Beta=Beta1)

    result


}
\keyword{ package }
