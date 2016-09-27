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

depthThreshold <- 3
offspringThreshold <- 200



#####

# find annotation
genes <- unique(as.character(annot[,3]))
#conceptsMat<-matrix(nrow=length(genes), ncol=20)
conceptsMat<-matrix(nrow=5000, ncol=20)
conceptsList<-list()
ctr<-1
# go through all genes and find terms from 
for(i in 1:50) {
#for(i in 1:length(genes)) {
  annotTerms<-annot[annot[,3]==genes[i] & annot[,9]==ontology & annot[,7]!="IEA",]
  goTerms<-unique(as.character(annotTerms[,5]))
  if (length(goTerms)<2) 
    next
  depths<-sapply(goTerms, findDepth)
  if (max(depths)<2) next;
  
  concepts<-vector()
  # create a square matrix of terms
  tmplca <- matrix(nrow = length(goTerms), ncol=length(goTerms))
  tmpdepth <- matrix(nrow = length(goTerms), ncol=length(goTerms))
  rownames(tmplca)<-goTerms
  colnames(tmplca)<-goTerms
  rownames(tmpdepth)<-goTerms
  colnames(tmpdepth)<-goTerms
  for(r in 1:nrow(tmplca)) {
      for(c in 1:ncol(tmplca)) {
        if (r==c) {
          # find lca 
          tmplca[r,c] <- NA
          # find depth of lca
          tmpdepth[r,c]<- -1
        }
        else {
          # find lca 
          tmplca[r,c] <- findLCA(goTerms[r], goTerms[c])
          # find depth of lca
          tmpdepth[r,c]<- findDepth(tmplca[r,c])
        }
      }
    }
  # loop through all terms to generalize concepts
  while(length(goTerms)>0){
    if(nrow(tmpdepth)>1) {
      if(max(tmpdepth) > depthThreshold) {
        # find max value in the matrix
        maxrow<-which(tmpdepth == max(tmpdepth), arr.ind = TRUE)[1,1]
        maxcol<-which(tmpdepth == max(tmpdepth), arr.ind = TRUE)[1,2]
        maxlca <- tmplca[maxrow, maxcol]
        concepts <- c(concepts, goTerms[maxrow], goTerms[maxcol])
        # remove max row and col
        browser()
        tmpdepth<-tmpdepth[-c(maxrow, maxcol),]
        tmpdepth<-tmpdepth[,-c(maxrow, maxcol)]
        # remove max row and col
        tmplca<-tmplca[-c(maxrow, maxcol), ]
        tmplca<-tmplca[ ,-c(maxrow, maxcol)]
        # check if lca exists in goTerms, if not add
        goTerms <- goTerms[-c(maxrow, maxcol)]
        if (!tmplca[maxrow,maxcol] %in% goTerms) {
          goTerms <- c(goTerms, maxlca)
        }
        MaxTermLCA <-sapply(goTerms, FUN=function(x) findLCA(x,maxlca))
        MaxTermLCADepth <- sapply(MaxTermLCA, findDepth)
        #MaxTermLCA <- c(MaxTermLCA, NULL)
        MaxTermLCADepth <- c(MaxTermLCADepth, -1)
        tmplca<-rbind(tmplca,MaxTermLCA)
        tmplca<-cbind(tmplca,c(MaxTermLCA,NULL))
        tmpdepth<-rbind(tmpdepth,MaxTermLCADepth)
        tmpdepth<-cbind(tmpdepth,c(MaxTermLCADepth,NULL))
      }
      # else no more term with depth greater than threshold
      else {
        concepts<-c(concepts, goTerms)
        goTerms<-NULL
      }
    }
    else {
      concepts<-c(concepts, goTerms)
      goTerms<-NULL
    }
  } #while(length(goTerms)>0)
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
  rm(concepts)
  rm(tmplca)
  rm(tmpdepth)
  ctr<-ctr+1
}
write.table(conceptsMat,file=paste("../concepts_",ont,"_thres_",depthThreshold,".txt",sep=""))
# perform association analysis on list
trans<-as(conceptsList,"transactions")
saveRDS(trans, file=paste("../trans_",ont,"_thres_",depthThreshold,".RDS",sep=""))
summary(trans)

