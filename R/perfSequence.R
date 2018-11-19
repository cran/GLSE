perfSequence <-
function(X){
  try(if( is.matrix(X)==FALSE)stop("X is not a matrix"))
  column<-dim(X)[[1]]
  dimnames(X) = list(c(seq(1:column)),c(seq(1:column)))
  diag(X) <- 0
  X[X!=0]=1
  Decomposable_test <- mcs(X)
  try(if( length(Decomposable_test) ==0)stop("X is not a decomposable graph"))
  chest.rip <- rip(X)
  return("Perfect.sequence"=chest.rip)
}
