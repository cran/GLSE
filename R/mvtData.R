mvtData <-
function(sigma=sigma, n,mean, df=1,delta = rep(0, nrow(sigma)),type = c("shifted", "Kshirsagar"),sd,Beta=NULL){
  N<-n
  column<-dim(sigma)[[1]]
  if(is.null(Beta)==TRUE){
    X <- rmvt(n=N, sigma=sigma, df=df)
    X<-data.frame(X)
    dimnames(X)[2]<-list(c(seq(1:column)))
    return1<-cbind(X=X)}

  else{
    try(if( is.vector(Beta)==FALSE)stop("Beta is not a vector"))
    try(if((dim(sigma)[[1]])!=length(Beta))stop(" Beta and sigma have non-conforming size"))
    X <- rmvt(n=N, sigma=sigma, df=df)
    bB<-Beta
    error<-rnorm(N, mean = mean, sd =sd)
    y<-X%*%bB+error
    colnames(y) <- "Y"
    X<-data.frame(X)
    dimnames(X)[2]<-list(c(seq(1:column)))
    return1<-cbind(X=X,y)}
  return("Data"=return1)
}
