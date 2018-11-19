mvnData <-
function(sigma, n,mean, method= c("eigen", "svd", "chol"),sd=1,Beta=NULL){



  if(is.null(Beta)==TRUE){
    column<-dim(sigma)[[1]]
    N<-n
    mu <-  mean
    X <- rmvnorm(N, mean=mu, sigma=sigma,method=method)
    X<-as.data.frame(X)
    dimnames(X)[2]<-list(c(seq(1:column)))
    return1<-cbind(X=X)}
  else{

    try(if( is.vector(Beta)==FALSE)stop("Beta is not a vector"))
    try(if((dim(sigma)[[1]])!=length(Beta))stop(" Beta and sigma have non-conforming size"))
    column<-dim(sigma)[[1]]
    N<-n
    mu <-  mean
    X <- rmvnorm(N, mean=mu, sigma=sigma,method=method)
    bB<-Beta
    error<-rnorm(N, mean = mean, sd =sd)
    y<-X%*%bB+error
    colnames(y) <- "Y"
    X<-as.data.frame(X)
    dimnames(X)[2]<-list(c(seq(1:column)))
    return1<-cbind(X=X,y)}

  return("Data"=return1)
}
