library(RBGL)
library(GO.db)
library(GOSim)
library(GOSemSim)
library(GOstats)
library(RBGL)
library(arules)

source("helperFunctions.R")
genes <-read.table("../genes.txt")
annot<-read.csv(file ="../gene_association/gene_association.sgd.csv",header = F)

ontology<-"P";
setOntology("BP", loadIC=FALSE)

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

depthThreshold <- Inf
offspringThreshold <- 200



#####

# find annotation
genes <- unique(as.character(annot[,3]))
#conceptsMat<-matrix(nrow=length(genes), ncol=20)
conceptsMat<-matrix(nrow=5000, ncol=60)
conceptsList<-list()
ctr<-1
# go through all genes and find terms from 
for(i in 1:length(genes)) {
#for(i in 1:500) {
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
  lstVector<-vector()
  for(j in 1:length(concepts)) {
    conceptsMat[ctr,j+1]<-concepts[j]
    #conceptsList[[ctr]][j]<-concepts[j]
    lstVector<-c(lstVector,concepts[j])
  }
  conceptsList<-c(conceptsList,list(lstVector))
  names(conceptsList)[ctr]<-genes[i]
  rm(lstVector)
  ctr<-ctr+1
}
write.table(conceptsMat,file=paste("../concepts_",ont,"_thres_",depthThreshold,".txt",sep=""))
# perform association analysis on list
trans<-as(conceptsList,"transactions")
saveRDS(trans, file=paste("../trans_",ont,"_thres_",depthThreshold,".RDS",sep=""))
summary(trans)
## generating rules
rules<-apriori(trans,parameter = list(supp = 0.001, conf = 0.8))
inspect(rules)
rules_high_lift<-sort(rules,by="lift")

