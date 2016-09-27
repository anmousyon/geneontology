library(RBGL)
library(GO.db)
library(GOSim)
library(GOSemSim)
library(GOstats)
library(RBGL)

genes <-read.table("../genes.txt")
annot<-read.csv(file ="../gene_association/gene_association.sgd.csv",header = F)

ontology<-"F";
setOntology("MF", loadIC=FALSE)

if (ontology=='P') {
  rootTerm<-"GO:0008150";
  ontoenv <- GOBPANCESTOR;
  parents <- GOBPPARENTS
  offspring <- as.list(GOBPOFFSPRING);
  ont ="BP"
} else if (ontology=="F") {
  rootTerm<-"GO:0003674"
  ontoenv <- GOMFANCESTOR
  parents <- GOMFPARENTS
  offspring <- as.list(GOMFOFFSPRING)
  ont = "MF"
}  else if (ontology=="C") {
  rootTerm<-"GO:0005575"
  ontoenv <- GOCCANCESTOR
  parents <- GOCCPARENTS
  offspring <- as.list(GOCCOFFSPRING)
  ont = "CC"
}

ancestors<-ontoenv

depthThreshold <- 3
offspringThreshold <- 200



#####

# find annotation
genes <- unique(as.character(annot[,3]))
#conceptsMat<-matrix(nrow=length(genes), ncol=20)
conceptsMat<-matrix(nrow=5000, ncol=20)
ctr<-1
# go through all genes and find terms from 
#for(i in 1:length(genes)) {
for(i in 1:500) {
  annotTerms<-annot[annot[,3]==genes[i] & annot[,9]==ontology & annot[,7]!="IEA",]
  goTerms<-unique(as.character(annotTerms[,5]))
  if (length(goTerms)<2) 
	  next
  depths<-sapply(goTerms, findDepth)
  if (max(depths)<2) next;
	  
  concepts<-vector()
  # loop through all terms to generalize concepts
  while(length(goTerms)>0){
    # always start with first term and then remove when done
    t<-goTerms[1]
    if (length(goTerms)>1) {
      otherTerms <- goTerms[-which(goTerms==t)]
      lca<-as.character(sapply(otherTerms, FUN=function(k) findLCA(t,k)))
      lcadepth<-sapply(lca, findDepth)
      lcaOffspringCount<-sapply(lca, findOffspringCount)
      # find lca with max depth
      if (max(lcadepth) > depthThreshold) {
        maxIndex <- which.max(lcadepth)
        # remove the first and the closest term
        goTerms<-goTerms[-c(1,maxIndex+1)]
        # replace with their lca if not already there
        if (!lca[maxIndex] %in% goTerms)
          goTerms<-c(lca[maxIndex],goTerms)
      }
      else {
        concepts<-c(concepts,goTerms[1])
        goTerms<-goTerms[-1]
      }
    }
    else { # only one term left 
      concepts<-c(concepts, goTerms[1])
      goTerms<-goTerms[-1]
    }
  }
  # write concepts vector in a file
  conceptsMat[ctr,1]<-genes[i]
  for(j in 1:length(concepts))
    conceptsMat[ctr,j+1]<-concepts[j]
  ctr<-ctr+1
  rm(concepts)
}
write.table(conceptsMat,file=paste("../concepts_",ont,".txt",sep=""))



 
findLCA<-function(term1, term2) {
  tryCatch ({
  if (length(get(term1,ontoenv))>0 && length(get(term2,ontoenv))>0) {
    g1 <- oneGOGraph(term1, parents)
    sg1 <- subGraph(c(get(term1, ontoenv), term1), g1)
    g2 <- oneGOGraph(term2, parents)
    sg2 <- subGraph(c(get(term2, ontoenv), term2), g2)
    g3 <- join(sg1, sg2)
    anc1<-union(ancestors[[term1]],term1)
    anc2<-union(ancestors[[term2]],term2)
    commonAncestors <- intersect(anc1, anc2)
    commonAncestors<-commonAncestors[-which(commonAncestors=="all")]
    dist1 = dijkstra.sp(sg1,term1)$distances[commonAncestors]
    dist2 = dijkstra.sp(sg2,term2)$distances[commonAncestors]
    dist <- dist1 + dist2
    if (length(dist)>0) 
      lca <- names(dist[which(dist==min(dist))])
    if (length(lca) > 1) {
      # find lowest lca i.e. with greatest depth
      depths <- sapply(lca, findDepth)
      maxind <- which(depths==max(depths))[1]
      lca[maxind]
    } else {
      lca
    }
  }
  else
    rootTerm
  }, error=function(e) {
    rootTerm
  })
}

findDepth<-function(term) {
  tryCatch ({
  g1 <- oneGOGraph(term, parents)
  sg <- subGraph(c(get(term, ontoenv), term), g1)
  depth <-  dijkstra.sp(sg, term)$distances[rootTerm]
  as.integer(depth)
  }, error=function(e){
      -1
  })
}

findOffspringCount<-function(term) {
  count<- offspring[term]
  if(!is.na(count))
    length(unlist(count))
  else
    0
}