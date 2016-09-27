findLCA<-function(term1, term2) {
  tryCatch ({
  ##if both of the genes exist?
  if (length(get(term1,ontoenv))>0 && length(get(term2,ontoenv))>0) {
    ##get gene from graph
    g1 <- oneGOGraph(term1, parents)
    sg1 <- subGraph(c(get(term1, ontoenv), term1), g1)
    ##get next gene from graph
    g2 <- oneGOGraph(term2, parents)
    sg2 <- subGraph(c(get(term2, ontoenv), term2), g2)
    ##create a pseudo-gene using subgraphs of last two?
    g3 <- join(sg1, sg2)
    ##find all ancestors of gene one
    anc1<-union(ancestors[[term1]],term1)
    ##fins all ancestors of gene two
    anc2<-union(ancestors[[term2]],term2)
    ##find common ancestors for each gene
    commonAncestors <- intersect(anc1, anc2)
    ##subtract all ancestors from list which are chared by all
    commonAncestors<-commonAncestors[-which(commonAncestors=="all")]
    ##get the depth of the gene in the subgraph
    dist1 = dijkstra.sp(sg1,term1)$distances[commonAncestors]
    ##get the depth of the gene in the subgraph
    dist2 = dijkstra.sp(sg2,term2)$distances[commonAncestors]
    ##combine their distances
    dist <- dist1 + dist2
    ##if they are not both the root
    if (length(dist)>0){ 
      ##get the least common ancestors
      lca <- names(dist[which(dist==min(dist))])
    }
    ##this could cause a problem if they are both root state since lca isnt defined it they are
    ##if they have more than one shared ancestor
    if (length(lca) > 1) {
      # find lowest lca i.e. with greatest depth
      depths <- sapply(lca, findDepth)
      maxind <- which(depths==max(depths))[1]
      lca[maxind]
    }
    ##otherwise use the only one we have (what if they dont have one from an error?)
    else {
      lca
    }
  }
  ##otherwise use the root as the lca (guarranteed to be an ancestor of both by nature of graph)
  else
    rootTerm
  }, error=function(e) {
    rootTerm
  })
}

findDepth<-function(term) {
  tryCatch ({
  ##get gene from graph
  g1 <- oneGOGraph(term, parents)
  sg <- subGraph(c(get(term, ontoenv), term), g1)
  ##find depth of the gene in the graph using dijkstras alg
  depth <-  dijkstra.sp(sg, term)$distances[rootTerm]
  ##convert to integer
  as.integer(depth)
  }, error=function(e){
      -1
  })
}

findOffspringCount<-function(term) {
  ##count the number of children the gene has
  count<- offspring[term]
  if(!is.na(count))
    length(unlist(count))
  else
    0
}