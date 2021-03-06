##' @name Information Centrality Functions
##' @rdname info.centrality
##'
##' @title Extensions to iGraph for Information Centrality
##'
##' @description Functions to compute the information centrality of a vertex (node) and network respectively. Includes a network efficiency measure to compute as a metric for information centrality. Uses graphs functions as an extension of \code{\link[igraph]{igraph}}.
##'
##' @param graph An \code{\link[igraph]{igraph}} object. May be directed or weighted as long as a shortest path can be computed.
##' @param verbose Logical. Whether computing information centrality of each node prints to monitor progress of a potentially long run-time. Defaults to FALSE.
##' @param net Numeric. Efficiency of the Network without any nodes removed. Defaults to computing for Graph given as input, can be given as a numeric if computed in advance to save run time.
##' @keywords graph network igraph centrality
##' @import igraph
NULL

##' @rdname info.centrality
##' @examples
##'
##' #generate example graphs
##' library("igraph")
##' g1 <- make_ring(10)
##' g2 <- make_star(10)
##'
##' #show network paths
##' distances(g1)
##' shortest_paths(g1, 5)
##'
##' #compute efficiency of full graphs
##' network.efficiency(g1)
##' network.efficiency(g2)
##' @export
network.efficiency <- function(graph){
  if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
  dd <- 1/shortest.paths(graph)
  diag(dd) <- NA
  efficiency <- mean(dd, na.rm=T)
  #denom <- nrow(dd)*(ncol(dd)-1)
  #sum(dd, na.rm=T)/denom
  return(efficiency)
}

##' @rdname info.centrality
##' @examples
##'
##' #compute information centrality (relative efficency when removed) for each node
##' info.centrality.vertex(g1)
##' info.centrality.vertex(g2)
##' @export
info.centrality.vertex <- function(graph, net=NULL, verbose=F){
  if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
  if(is.null(net)) net <- network.efficiency(graph)
  if(is.numeric(net)==F){
    warning("Please ensure net is a scalar numeric")
    net <- network.efficiency(graph)
  }
  count <- c()
  for(i in 1:length(V(graph))){
    count <- c(count, (net-network.efficiency(delete.vertices(graph, i)))/net)
    if(verbose){
      print(paste("node",i,"current\ info\ score", count[i], collapse="\t"))
    }
  }
  return(count)
}
##' @rdname info.centrality
##' @examples
##'
##' #compute total information centrality for a network
##' info.centrality.network(g1)
##' info.centrality.network(g2)
##' @export
info.centrality.network <- function(graph, net=network.efficiency(graph), verbose=F) sum(info.centrality.vertex(graph))