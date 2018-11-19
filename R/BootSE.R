BootSE <-function(Graph,x,y,n.boot,sig.level= 0.05,centering=TRUE){

    column <- dim(x)[[2]]
    n <- dim(x)[[1]]
    sq <- matrix(0,column,column)
    x <- as.matrix(x)
    y <- as.matrix(y)
    dimnames(Graph) <- list(c(seq(1:column)),c(seq(1:column)))
    chest.rip <- rip(Graph)
    Cliques=chest.rip$cliques
    Separators=chest.rip$separators
    Separators<-Separators[lapply(Separators,length)>0]
    if(centering){
      x=scale(x, center = TRUE, scale = FALSE)
      y=scale(y, center = TRUE, scale = FALSE)
    }
    samplecovmatrix=t(x)%*%x

    sum_Cliq <- sq

    for(l in 1:length(Cliques)){
      C <-  as.numeric(c(do.call("cbind",Cliques[l])) )
      sum_Cliq[C,C] <- sum_Cliq[C,C]+ chol2inv(chol((samplecovmatrix[C,C])))
    }
    sum_Sep <- sq
    if (length(Separators)>0){
      for(p in 1:length(Separators)){
        S <-  as.numeric(c(do.call("cbind",Separators[p])))
        sum_Sep[S,S] <- sum_Sep[S,S]+chol2inv(chol(samplecovmatrix[S,S]))
      }
    }
    beta <- (sum_Cliq-sum_Sep)%*%t(x)%*% y

    B <- list()
    i <- 0
    tt <- 0
    repeat{
      i <- i+1
      print(i)
      ### draw the bootstrap sample and calculate thetahat.b
      index = 1:n
      bootindex = sample(index, n, replace=TRUE)
      samplex = x[bootindex,]
      sampley = y[bootindex]
      xx <- as.matrix(samplex)
      yy <- as.matrix(sampley)
      samplecovmatrix=t(xx)%*%xx

      sum_Cliq <- sq
      for(l in 1:length(Cliques)){
        C <-  as.numeric(c(do.call("cbind",Cliques[l])) )
        if(det(as.matrix(samplecovmatrix[C,C]))<0.0001) break
        else{
          sum_Cliq[C,C] <- sum_Cliq[C,C]+ chol2inv(chol((samplecovmatrix[C,C])))
        }
      }
      sum_Sep <- sq
      if (length(Separators)>0){
        for(p in 1:length(Separators)){
          S <-  as.numeric(c(do.call("cbind",Separators[p])))
          if(det(as.matrix(samplecovmatrix[S,S]))<0.0001) break
          sum_Sep[S,S] <- sum_Sep[S,S]+chol2inv(chol(samplecovmatrix[S,S]))
        }
      }


      B[[i]] <- (sum_Cliq-sum_Sep)%*%t(xx)%*% yy
      if(any(as.vector(B[[i]])==0))
        B[[i]] <- NULL
      B[sapply(B, is.null)] <- NULL
      B.clean <- B
      if (length(B.clean)==n.boot) break
    }

    store.matrix <- matrix(unlist(B),nrow=column,ncol=length(B))
    SE <-  apply(store.matrix,1,sd)
    dec <- c()
    for(i in 1:column){
      qq <- quantile(store.matrix[i,],probs=c((sig.level/2),(1-(sig.level/2))))
      dec[i] <- ifelse(qq[[1]]<0&qq[[2]]>0,"+","**")
    }
    return(data.frame("Beta"=round(beta,digits=4),"Sttandard.Error"=round(SE,digits=4),"Significance"=dec))
}
