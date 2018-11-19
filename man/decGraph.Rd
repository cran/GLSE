\name{decGraph}
\alias{decGraph}
\title{
  Generating a decomposable graph
}
\description{ The function helps in generating a decomposable graph based on the given number of vertices and edges in the graph. It supports both serial and parallel computations. However, the available computational resources should be considered before choosing one of the two options. The return will be a decompsable graph. 
}
\usage{
decGraph(X, edge, parallel = FALSE, clusters = NULL)
}
\arguments{
    \item{X}{
      The number of vertices required in the generated decomposable graph. It should be noted that the function will be slow if the number of vertices is more than 50.
    }
    \item{edge}{
      The number of edges required in the generated decomposable graph (the maximum number of edges \code{(v*(v-1)/2)-1}).
    }
    \item{parallel}{
     Logical. Whether or not to do parallel computation. If set to \code{TRUE}, then \code{clusters} must also be provided. The default is \code{FALSE}.
    }
    \item{clusters}{
      The number of child nodes that must be created in the cluster if the \code{parallel} option is set to \code{TRUE}.
    }
  }

  \value{
The output is
    \item{Perfect.sequence}{ The perfect sequence of cliques and separators for a given graph.}

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


###################################

# Example with serial processing.

###################################
 library(igraph)
library(gRbase)

# specified the number of vertices and edges

X=10
edge= 20

# Generate a decomposable graph

Graph <-decGraph (X,edge,parallel=FALSE,clusters=NULL)

###################################

# Example with parallel processing.

###################################
library(parallel)
# Create cluster

cl <- makeCluster(2)

# Send the gRbase package to all the nodes

clusterEvalQ(cl, library(gRbase))


# Generate a decomposable graph.

Graph <-decGraph (X,edge,parallel=TRUE,clusters=cl)

# Stop the cluster

stopCluster(cl)

Graph



  }
