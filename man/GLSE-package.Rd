\name{GLSE-package}
\alias{GLSE-package}
\docType{package}
\title{Graphical Least Square Estimation}
\description{
The Graphical Least Square Estimation 'GLSE' package consists of several functions. The first function (called GLSE()), which is the main function in this package, is for estimating linear regression via graphical models, especially when the number of observations is smaller than the number of variables. This function gives unbiased estimated regression models under some conditions on the predictors. Moreover, this function aims to deal with studies where all covariates need to be kept in the model (such as portfolio optimisation). The second function provides an estimation of the bootstrap standard errors of the 'GLSE' parameters. The third function generates data from multivariate normal distribution for regression purposes, and the fourth function generates data from multivariate-distribution for the same purpose. Generating an undirected decomposable graph is carried out by the fifth function and its decomposability is tested in the sixth function. The seventh function gives a perfect sequence of the cliques and separators of the given decomposable graph. Finally, the last function gives the plot of the decomposale graph. The package supports both serial and parallel computations for estimating the graph.The 'GLSE' package is an implementation of the method introduced in Aldahmani, S., and Dai, H. (2015) <doi:10.5539/ijsp.v4n3p61>.  

}
\details{
\tabular{ll}{
Package: \tab GLSE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2018-07-04\cr
License: \tab GPL-2 \cr
}

}
\author{Saeed Aldahmani 

Maintainer: Saeed Aldahmani <saldahmani@uaeu.ac.ae>

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
\keyword{ package }



