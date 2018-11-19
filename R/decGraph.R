decGraph <-function(X,edge,parallel=FALSE,clusters=NULL){
  
  mat<-matrix(0,X,X)

  try(if(((X*(X-1)/2)-1)<edge)stop("too many edges"))
  indplace0 <- function(mat) which(mat==0 & upper.tri(mat, diag = F) == TRUE ,arr.ind=T)

  ff <- function(xxx){
    ind <- xxx
    mat[ind[1],ind[2]] <- 1
    ind11 <- lower.tri(mat)
    mat[ind11] <- t(mat)[ind11]
    dimnames(mat) <- list(c(seq(1:X)),c(seq(1:X)))
    Decomposable_test <- mcs(mat)
    if (length(Decomposable_test) >0){
      return("mat"=mat)}
    else{return("mat"=NA)}
  }
  if(parallel==FALSE){
    i=0
    matlist<-list()
    repeat{
      i=i+1
      if (i==1){
        ind1 <- indplace0(mat)
        ll22  <- split(ind1, rep(1:nrow(ind1), each = 1))
        matlist <- lapply(ll22,ff) }

      else{
        ind1 <- indplace0(mat)
        ll22  <- split(ind1, rep(1:nrow(ind1), each = 1))
        matlist <- lapply( ll22,ff)}

      if(any(is.na(matlist))==TRUE ){
        matlist <- matlist[-which(is.na(matlist))]
        matlist <- na.omit(matlist)}
      if (length(matlist)!=0){
        mat<-sample(matlist,1)[[1]]
        matt11<- mat}
      if(length(matlist)==0||i==edge) break
    }}
  else{
    i=0
    matlist<-list()
    cl <- clusters
    repeat{
      i=i+1
      if (i==1){
        ind1 <- indplace0(mat)
        ll22  <- split(ind1, rep(1:nrow(ind1), each = 1))
        matlist <- parLapply(cl,ll22,ff)  }

      else{
        ind1 <- indplace0(mat)
        ll22  <- split(ind1, rep(1:nrow(ind1), each = 1))
        matlist <- parLapply(cl,ll22,ff) }

      if(any(is.na(matlist))==TRUE ){
        matlist <- matlist[-which(is.na(matlist))]
        matlist <- na.omit(matlist)}
      if (length(matlist)!=0){
        mat<-sample(matlist,1)[[1]]
        matt11<- mat}
      if(length(matlist)==0||i==edge) break
    }}
  return("mat"= matt11)
}
