\name{testDeco}
\alias{testDeco}
\title{
  Test of undirected decomposable graph
}
\description{ The function examines if the graph is an undirected decomposable one.
}
\usage{
testDeco(X)
}
\arguments{
    \item{X}{
      A matrix of \code{v} x  \code{v} dimension.
    }

  }

  \value{
    The output consists of either
    \item{a decomposable graph}{ The given concentration matrix (graph) is decomposable.}
    \item{non-decomposable graph}{ The given concentration matrix (graph) is not decomposable.}
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
    \code{\link{decGraph}}
  }
    \examples{
library(igraph)
library(gRbase)

# The number of vertices and edges

X=10
edge= 20

# Generate a decomposable graph

Graph <-decGraph (X,edge,parallel=FALSE,clusters=NULL)

# Call testDec(.).

testDeco ( Graph)

  }
