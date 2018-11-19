testDeco <-
function(X){
  try(if( is.matrix(X)==FALSE)stop("X is not a matrix"))
  column<-dim(X)[[1]]
  dimnames(X) = list(c(seq(1:column)),c(seq(1:column)))
  diag(X) <- 0
  X[X!=0]=1
  Decomposable_test <- mcs(X)
  if(length(Decomposable_test) ==0){
    return1<-"Not decomposable graph"
  }
  else{ return1<-"Decomposable graph"}


  return(return1)
}
