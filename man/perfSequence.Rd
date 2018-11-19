\name{perfSequence}
\alias{perfSequence}
\title{
Perfect sequence of the graph
}
\description{ The function gives a perfect sequence of cliques and separators for a given decomposable concentration matrix.
}
\usage{
perfSequence(X)
}

\arguments{
  \item{X}{
A decomposable concentration matrix of dimension \code{v} x  \code{v}.
}

}

\value{
A data frame consisting of
 \item{Perfect.sequence}{ The perfect sequence of cliques and separators and parents for a given concentration matrix.}

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
 \code{\link{testDeco}}
}
\examples{
library(igraph)
library(gRbase)

# The number of vertices and edges

x=10
edge= 20

# Generate a decomposable graph

Graph <-decGraph (x,edge)

# Call perfSequence(.).


perfSequence ( Graph )



}

