library(RBGL)
library(GO.db)
library(GOSim)
library(GOSemSim)
library(GOstats)
library(RBGL)
library(arules)

source("helperFunctions.R")
genes <-read.table("../genes.txt")
annot<-read.csv(file ="../gene_association/human_association.txt",header = F)

ontology<-"F";
setOntology("MF", loadIC=FALSE)

##select the right configuration for the chosen ontology
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

##set ancestors
ancestors<-ontoenv

##initialize thresholds
depthThreshold <- 4
offspringThreshold <- 200



#####

# find annotation
##get fthird column of gene annotation and store in list (set since unique)
genes <- unique(as.character(annot[,3]))
#conceptsMat<-matrix(nrow=length(genes), ncol=20)
##create a matrix of size 5000*60
conceptsMat<-matrix(nrow=5000, ncol=60)
##create an empty list to store the concepts
conceptsList<-list()
##?
ctr<-1
# go through all genes and find terms from
##for all genes in list
#for(i in 1:length(genes)) {
for(i in 1:500) {
  ##get all the information for that gene, for specific ontology
  annotTerms<-annot[annot[,3]==genes[i] & annot[,9]==ontology & annot[,7]!="IEA",]
  ##get fifth column of annotation and store in list (set since unique)
  goTerms<-unique(as.character(annotTerms[,5]))
  ##if there are less than two terms then skip this gene
  if (length(goTerms)<2) 
    next
  ##store the depths of all the terms in a list
  depths<-sapply(goTerms, findDepth)
  ##if none of the terms have a depth greater than 2 then skip this gene
  if (max(depths)<2) next;
  
  ##create a vector to store the concepts
  concepts<-vector()
  # loop through all terms to generalize concepts
  while(length(goTerms)>0){
    # always start with first term and then remove when done
    t<-goTerms[1]
    ##if there are at least two terms left
    if (length(goTerms)>1) {
      ##otherterms is all terms that are not the same as the first one (t)?
      otherTerms <- goTerms[-which(goTerms==t)]
      ##find lca by using the findLCA helper function
      lca<-as.character(sapply(otherTerms, FUN=function(k) findLCA(t,k)))
      ##find depth using lcaDepth helper function
      lcadepth<-sapply(lca, findDepth)
      ##find offspring using lcaOffspringCount helper function
      lcaOffspringCount<-sapply(lca, findOffspringCount)
      # find lca with max depth
      ##if any of the terms have a depth greater than the threshold
      if (max(lcadepth) > depthThreshold) {
        ##store greatest term depth
        maxIndex <- which.max(lcadepth)
        # remove the first and the closest term
        ##remove first term using -1 and closest using -(maxindex+1)
        goTerms<-goTerms[-c(1,maxIndex+1)]
        # replace with their lca if not already there
        ##if the term with the max depth isnt in goTerms
        if (!lca[maxIndex] %in% goTerms)
          ##add it to the goTerms list?
          goTerms<-c(lca[maxIndex],goTerms)
      }
      ##if no terms have a depth greater than the threshold
      else {
        ##store the concept of the first term
        concepts<-c(concepts,goTerms[1])
        ##remove the first element from the list
        goTerms<-goTerms[-1]
      }
    }
    ##if there is only one term
    else { # only one term left 
      ##store the concept of the first term
      concepts<-c(concepts, goTerms[1])
      ##remove the first element from the list
      goTerms<-goTerms[-1]
    }
  }
  # write concepts vector in a file
  ##store gene i in row=ctr column=1
  conceptsMat[ctr,1]<-genes[i]
  ##create a new vector
  lstVector<-vector()
  ##for every concepyt
  for(j in 1:length(concepts)) {
    ##insert the concept into matrix at row=ctr column = j+1
    ##format will be [gene, c1, c2, c3, c4, ...]
    conceptsMat[ctr,j+1]<-concepts[j]
    #conceptsList[[ctr]][j]<-concepts[j]
    ##append the concept to lstVector
    lstVector<-c(lstVector,concepts[j])
  }
  ##append the list of concepts to the concept list
  ##will be list of lists?
  conceptsList<-c(conceptsList,list(lstVector))
  names(conceptsList)[ctr]<-genes[i]
  rm(lstVector)
  ##increment the row number?
  ctr<-ctr+1
}
write.table(conceptsMat,file=paste("../concepts_",ont,"_thres_",depthThreshold,".txt",sep=""))
# perform association analysis on list
##get all the transactions for the concepts
trans<-as(conceptsList,"transactions")
saveRDS(trans, file=paste("../trans_",ont,"_thres_",depthThreshold,".RDS",sep=""))
summary(trans)
## generating rules
rules<-apriori(trans,parameter = list(supp = 0.001, conf = 0.8))
inspect(rules)
##sort the rules by their lift value
rules_high_lift<-sort(rules,by="lift")
write(rules_high_lift, file=paste("../rules_",ont,"_thresh_",depthThreshold,"_sort_lift.txt", sep=""), sep="\t", col.names=NA)
