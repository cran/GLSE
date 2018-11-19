\name{BootSE}
\alias{BootSE}
\title{
Bootstrap  error estimation	and	significance	for	the	GLSE 	estimates.
}
\description{
The function provides bootstrap error estimation and significance for the estimates obtained via the
GLSE().
}
\usage{
BootSE(Graph, x, y, n.boot, sig.level = 0.05, centering = TRUE)
}

\arguments{
  \item{Graph}{
An object of class \code{GLSE}.
}
  \item{x}{
A matrix of \code{n} x  \code{v} dimension consisting of the predictors values.
}
  \item{y}{
A vector of length \code{n} of the response variable.
}
  \item{n.boot}{
Number of bootstrap samples.
}
  \item{sig.level}{
Significance level. Must be pre-specified.
}
 \item{centering}{
Logical. If equal to \code{TRUE} then all the variables in \code{x} and the response \code{y} are centered. The default is \code{TRUE}.
}
}

\value{
A data frame consisting of
 \item{Beta }{ The estimated regression coefficients based on the provided graph.}
 \item{Sttandard.Error }{ Bootstrap estimates of the standard errors.}
\item{Significance }{ Significance of the estimates based on the specified quantile points. \code{"**"} mean significance at the specified level while \code{"+"} means no significance.}
}
\references{
Alan, G. Frank, B. (2009), Computation of Multivariate Normal and t Probabilities.
  Lecture Notes in Statistics, Vol. 195., Springer-Verlag, Heidelberg. ISBN 978-3-642-01688-2.
  
   Alan, G. Frank, B. Tetsuhisa, M. Xuefe, M. Friedrich, L. Fabian, S. and Torsten,
  H. (2019). mvtnorm: Multivariate Normal and t Distributions. R package version 1.0-7.
  URL http://CRAN.R-project.org/package=mvtnorm.
  
    Aldahmani, S. and Dai, H. (2015). Unbiased Estimation for Linear Regression When n< v. International Journal of Statistics and Probability, 4(3), p61.

 Csardi, G., and Nepusz, T. (2006). The igraph software package for complex network research. InterJournal, Complex Systems, 1695(5), 1-9.
 
    Dethlefsen, C., and H?jsgaard, S. (2005). A common platform for graphical models in R: The gRbase package. Journal of Statistical Software, 14(17), 1-12.}
    
  
\author{
Saeed Aldahmani 
}

\seealso{
 \code{\link{GLSE}}
}

\examples{
library(gRbase)
#load the data
data(qsbralks)
x <- as.matrix(qsbralks[,1:15])
y <- as.matrix(qsbralks[,22])

# Call GLSE(.) .


Graph <- GLSE(x,y,parallel=FALSE)$Graph


# Call BootSE(.) .

results <- BootSE(Graph,x,y,n.boot=150,sig.level =.3,centering = TRUE)

results

}

