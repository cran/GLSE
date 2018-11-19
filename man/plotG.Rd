\name{plotG}
\alias{plotG}
\title{
  Plotting an undirected decomposable graph
}
\description{ The function plots a given undirected decomposable graph.}

\usage{
 plotG(X, vertex_size = 15, vertex_color = "green",
 vertex_frame_color = "green", vertex_label_cex = 1, edge_width = 3)
}
 \arguments{
    \item{X}{
      A concentration matrix of \code{v} x  \code{v} dimension.
    }
    \item{vertex_size}{
      The size of the vertex, a numeric scalar or vector.
    }
    \item{vertex_color}{
      A vector of the node color which should specify a color name or number.
    }
       \item{vertex_frame_color}{
      A vector of the node frame color which should specify a color name or number.
    }
     \item{vertex_label_cex}{
      A numeric vector of node label size. }
    \item{edge_width}{
     Edge width.
    }
  }

  \value{
    The output is
    \item{an undirected graph}{The plotted graph of the given matrix.}

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

v=10
edge= 20

# Generate a decomposable graph

Graph <-decGraph (v,edge,parallel=FALSE,clusters=NULL)


# Call plotG(.).

plotG(Graph,vertex_size=15,vertex_color="green",
vertex_frame_color="green",vertex_label_cex=1, edge_width=3)


  }
