GLSE <-function(x,y,centering=TRUE,lambda=0,parallel=FALSE,clusters=NULL,c.size=NULL){

  column <- dim(x)[[2]]
  n <- dim(x)[[1]]
  if(is.null(c.size)){
    rr <- n
  }
  else{try(if( c.size > n) stop("c.size > n"))
    rr <- c.size}
  n.edge <- column*(column-1)/2
  x <- as.matrix(x)
  y <- as.matrix(y)
  #
  if(centering){
    x=scale(x, center = TRUE, scale = FALSE)
    y=scale(y, center = TRUE, scale = FALSE)
  }
  samplecovmatrix=t(x)%*%x
  cl <- clusters
  sq <- matrix(0,column,column)

  indplace0 <- function(mat) which(mat==0 & upper.tri(mat, diag = FALSE) == TRUE ,arr.ind=TRUE)
  if(parallel==FALSE){
    f1 <- function(xx,k){
      edges <- k
      mat  <- xx
      ind1 <- indplace0(mat)
      ll1  <-  split(ind1, rep(1:nrow(ind1), each = 1))
      fun1 <- function(newmat){
        dimnames(newmat) <- list(c(seq(1:column)),c(seq(1:column)))
        chest.rip <- rip(newmat)
        if (!is.null(chest.rip)){
          Cliques=chest.rip$cliques
          if(rr >max(sapply(Cliques,length))){
            sum_Cliq <- sq
            for(l in 1:length(Cliques)){
              C <-  as.numeric(c(do.call("cbind",Cliques[l])) )
              sum_Cliq[C,C] <- sum_Cliq[C,C]+ chol2inv(chol((samplecovmatrix[C,C])))
            }
            sum_Sep <- sq
            Separators=chest.rip$separators
            Separators<-Separators[lapply(Separators,length)>0]
            if (length(Separators)>0){
              for(p in 1:length(Separators)){
                S <-  as.numeric(c(do.call("cbind",Separators[p])))
                sum_Sep[S,S] <- sum_Sep[S,S]+chol2inv(chol(samplecovmatrix[S,S]))
              }
            }
            SSE1 <- (sum((y- (x%*%((sum_Cliq-sum_Sep)%*%t(x)%*% y)))^2)) + lambda*edges
          }
        }
        if (is.null(chest.rip)| !exists("SSE1"))
          return(NA)
        else
          return(SSE1)
      }
      ################################################
      ff <- function(xxx){
        ind <- xxx
        mat[ind[1],ind[2]] <- 1
        ind11 <- lower.tri(mat)
        mat[ind11] <- t(mat)[ind11]
        S <- fun1(mat)
        return(list("mat"=mat,"error"=S))
      }
      matlist <- lapply(ll1,ff)
      return(matlist)
    }
  }
  else{
    f1 <- function(xx,k){

      edges <- k
      mat  <- xx
      ind1 <- indplace0(mat)
      ll1  <-  split(ind1, rep(1:nrow(ind1), each = 1))
      fun1 <- function(newmat){
        dimnames(newmat) <- list(c(seq(1:column)),c(seq(1:column)))
        chest.rip <- rip(newmat)
        if (!is.null(chest.rip)){
          Cliques=chest.rip$cliques
          if(n >max(sapply(Cliques,length))){
            sum_Cliq <- sq
            for(l in 1:length(Cliques)){
              C <-  as.numeric(c(do.call("cbind",Cliques[l])) )
              sum_Cliq[C,C] <- sum_Cliq[C,C]+ chol2inv(chol((samplecovmatrix[C,C])))
            }
            sum_Sep <- sq
            Separators=chest.rip$separators
            Separators<-Separators[lapply(Separators,length)>0]
            if (length(Separators)>0){
              for(p in 1:length(Separators)){
                S <-  as.numeric(c(do.call("cbind",Separators[p])))
                sum_Sep[S,S] <- sum_Sep[S,S]+chol2inv(chol(samplecovmatrix[S,S]))
              }
            }
            SSE1 <- (sum((y- (x%*%((sum_Cliq-sum_Sep)%*%t(x)%*% y)))^2)) + lambda*edges
          }
        }
        if (is.null(chest.rip)| !exists("SSE1"))
          return(NA)
        else
          return(SSE1)
      }
      ################################################
      ff <- function(xxx){
        ind <- xxx
        mat[ind[1],ind[2]] <- 1
        ind11 <- lower.tri(mat)
        mat[ind11] <- t(mat)[ind11]
        S <- fun1(mat)
        return(list("mat"=mat,"error"=S))
      }
      matlist <- parLapply(cl,ll1,ff)
      return(matlist)
    }
  }
  k=0
  M <- list()
  tot=0
  ssef <- c()
  repeat{
    k <- k+1
    


    if(k==1){
      mat <- sq
      result <- f1(mat,k)
      sse <- sapply(result,"[[","error")
      matlist <- sapply(result,"[","mat")
    }
    else{
      result <- f1(mat,k)
      sse <- sapply(result,"[[","error")
      matlist <- sapply(result,"[","mat")
    }
    if(any(!is.na(sse))==FALSE) break
    if(k==(column-1)) break
    else{
      if(any(is.na(sse))==TRUE){
        matlist <- matlist[-which(is.na(sse))]
        sse  <- na.omit(sse)
      }

      mat <- matlist[[which.min(sse)]]
      M[[k]] <- mat
      ssef[k] <- min(sse)
    }
  }
  bm <- M[[which.min(ssef)]]
  dimnames(bm) <- list(c(seq(1:column)),c(seq(1:column)))
  chest.rip <- rip(as.matrix(bm))

  Cliques=chest.rip$cliques
  sum_Cliq<- sq
  for(l in 1:length(Cliques)){
    C <-  as.numeric(c(do.call("cbind",Cliques[l])) )

    sum_Cliq[C,C] <- sum_Cliq[C,C]+ chol2inv(chol((samplecovmatrix[C,C])))

  }
  sum_Sep<- sq
  Separators=chest.rip$separators
  Separators<-Separators[lapply(Separators,length)>0]
  if (length(Separators)>0){
    for(k in 1:length(Separators)){
      S <-  as.numeric(c(do.call("cbind",Separators[k])))
      sum_Sep[S,S] <- sum_Sep[S,S]+chol2inv(chol((samplecovmatrix[S,S])))
    }
  }
  g <- sum_Cliq-sum_Sep
  BO <- (g)%*%t(x)%*% y
  SSE <- sum((y- (x%*%BO))^2)
  Error <- sqrt(diag(((g%*%t(x))%*% t(g%*%t(x)))*(SSE/(n)))+diag(((g%*%t(x)%*%x%*%BO)-BO)%*%t(((g%*%t(x)%*%x%*%BO)-BO))))

  return(list("Graph"=bm,"Beta"=BO,"SSE"=SSE,"Error"=Error))
}