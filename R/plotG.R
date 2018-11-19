plotG <-
function(X,vertex_size=15,vertex_color="green",vertex_frame_color="green",vertex_label_cex=1, edge_width=3){
 
   try(if( is.matrix(X)==FALSE)stop("X is not a matrix"))
  column<-dim(X)[[1]]
  diag(X) <- 0
  X[X!=0]=1
  adjacency.plot <- graph.adjacency(X, mode='undirected')
  pl<-plot(adjacency.plot,vertex.size=vertex_size,vertex.color=vertex_color,vertex.frame.color=vertex_frame_color,vertex.label.cex=vertex_label_cex, edge.width=edge_width)
}
