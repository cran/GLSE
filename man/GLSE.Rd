\name{GLSE}
\alias{GLSE}
\title{
Graphical least square estimation
}
\description{
The function estimates linear regression via graphical models in an unbiased way. This function aims to deal with studies where all covariates need to be kept in the model (such as portfolio optimisation), when the sample size is much less than covariate numbers. It supports both serial and parallel computations. However, the available computational resources should be considered before choosing one of the two options. For any given data, an estimated graph, regression parameters and errors are returned.
}
\usage{
GLSE(x, y, centering = TRUE, lambda = 0, parallel = FALSE, clusters = NULL, c.size = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix of \code{n} x  \code{v} dimension consisting of the predictors' values.
}
  \item{y}{
A response vector of length  \code{n}.
}
  \item{centering}{
If equal to \code{TRUE}, then all the variables in \code{x} and the response \code{y} are centered. The default is \code{TRUE}.
}
  \item{lambda}{
The amount of regularisation for shrinking the graph. The default is \code{0}.
}
  \item{parallel}{Logical.Whether or not to do parallel computation. If set to \code{TRUE}, then \code{clusters} must also be provided. The default is \code{FALSE}.
}
  \item{clusters}{
 The number of child nodes that must be created in the cluster if the \code{parallel} option is set to \code{TRUE}.
}
 

 \item{c.size}{
The allowed size of the clique (should not be greater than n). The default is \code{n}.
}
}
\details{
The parallel options also work for Windows. Care should be taken in creating the cluster as all types of clusters do not work for Windows, e.g. \code{FORK}.
}
\value{
A list of three values:
 \item{Graph }{ The estimated graph for the given data.}
  \item{Beta }{ The vector of estimated regression parameters based on the estimated graph.}
\item{SSE }{ Sum of squares of errors for the parameters.}
\item{SE }{ Estimated standard errors for the parameters.}

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
\author{
Saeed Aldahmani 
}
\seealso{
 \code{\link{BootSE}}
 }
\examples{

###################################

# Example with serial processing.

###################################
library(gRbase)
library(igraph)
#load the data
data(qsbralks)
x <- as.matrix(qsbralks[,1:15])
y <- as.matrix(qsbralks[,22])

# Get the graph

# Call GLSE(.) .

results <- GLSE(x,y,parallel=FALSE)


# The returned objects are

results$Graph
results$Beta
results$SSE
results$Error

# Call plotG(.).

 plotG(results$Graph)
###################################

# Example with parallel processing.

###################################
library(gRbase)
library(igraph)
library(parallel)
#load the data
data(qsbralks)
x <- qsbralks[,1:15]
y <- qsbralks[,22]

# Create cluster

cl <- makeCluster(2)

# Send the gRbase package to all the nodes

clusterEvalQ(cl, library(gRbase))


# Call GLSE(.) .

results <- GLSE(x,y,parallel=TRUE,clusters=cl)

# Stop the cluster

stopCluster(cl)

# The returned objects are

results$Graph;
results$Beta;
results$SSE;
results$Error;

# Call plotG(.) .

plotG(results$Graph)


}

\keyword{ package }
