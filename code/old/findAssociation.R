library(RBGL)
library(GO.db)
library(GOSim)
library(GOSemSim)
library(GOstats)
library(RBGL)

genes <-read.table("../genes.txt")
annot<-read.csv(file ="../gene_association/gene_association.sgd.csv",header = F)
genes_input<-read.table("../experiments/bma_approach_PhaseII/dataset1//BP//no_IEA//1//similarityP_final.txt")
########### parameters ######################################
NUM_TRIALS<-4955;
start<- 0;

OffspringTuningFactor <- 0.6;
ontology<-"C";
setOntology("CC", loadIC=FALSE)

if (ontology=='P') {
  rootTerm<-"GO:0008150";
  ontoenv <- GOBPANCESTOR;
  offspring <- as.list(GOBPOFFSPRING);
  ont ="BP"
} else if (ontology=="F") {
  rootTerm<-"GO:0003674"
  ontoenv <- GOMFANCESTOR
  offspring <- as.list(GOMFOFFSPRING)
  ont = "MF"
}  else if (ontology=="C") {
  rootTerm<-"GO:0005575"
  ontoenv <- GOCCANCESTOR
  offspring <- as.list(GOCCOFFSPRING)
  ont = "CC"
}

ancestors<-ontoenv
##############################################################

# find lca of group of terms


maxOffspring <- length(as.list(offspring)[rootTerm][[1]])
similarity<-matrix(nrow=NUM_TRIALS, ncol=15)
if (file.exists(paste("../similarity_cache_",ontology,".txt"))) {
  similarity_cache <- read.table(paste("similarity_cache_"),ontology,".txt")
} else {
similarity_cache<-matrix(nrow=0, ncol=4)
colnames(similarity_cache)<-c("term1","term2","ontology","sim")
}

for(ctr in 1:NUM_TRIALS) {
  gene1<-as.vector(genes[start+ctr,2]);
  gene2<-as.vector(genes[start+ctr,5]);
	#gene1<- as.vector(genes_input[ctr,1]);
	#gene2<- as.vector(genes_input[ctr,2]);
  #geneid1<-as.vector(genes[start+ctr,3])
  #geneid2<-as.vector(genes[start+ctr,6])
  #lookup annotations
  annot1<-annot[annot[,3]==gene1 & annot[,9]==ontology, ]  #& annot[,7]!="IEA",]
  goTerms1<-as.character(annot1[,5])
  annot2<-annot[annot[,3]==gene2 & annot[,9]==ontology, ] #& annot[,7]!="IEA",]
  goTerms2<-as.character(annot2[,5])
  if (length(goTerms1)==0 || length(goTerms2)==0)
    next;
  
  tempMatrixOffspring <-matrix(nrow=length(goTerms1), ncol=length(goTerms2))
  
  # find similarity between goTerms1 and goTerm2 
  for (i in 1:length(goTerms1)) {
    term1 <- goTerms1[i]
    g1 <- oneGOGraph(term1, ontoenv)
    for(j  in 1:length(goTerms2)) {
      term2 <- goTerms2[j]
      # check caching table to see if similarity exists
      if (length(similarity_cache[similarity_cache[,1]==term1 & similarity_cache[,2]==term2,]) >0)
      {
          simOffspring = similarity_cache[similarity_cache[,1]==term1 & similarity_cache[,2]==term2,4]
      }
      else {
        # find closest term from goTerms2
        g2 <- oneGOGraph(term2, ontoenv)
        g3 <- join(g1, g2)
        anc1<-union(ancestors[[term1]],term1)
        anc2<-union(ancestors[[term2]],term2)
        commonAncestors <- intersect(anc1, anc2)
        commonAncestors<-commonAncestors[-which(commonAncestors=="all")]
        dist1 = dijkstra.sp(g1,term1)$distances[commonAncestors]
        dist2 = dijkstra.sp(g2,term2)$distances[commonAncestors]
        dist <- dist1 + dist2
        # find min value of dist and how many terms have it
        if (length(dist)>0) {
          lcas <- names(dist[which(dist==min(dist))])
          if (length(lcas) > 1)
          {
            lcas <- names(sort(sapply(lcas, FUN=function(y) length(offspring[y][[1]]))))[1]
          }
        }
        else
          lcas <-NA
        
        minSubsumerIC <- getMinimumSubsumer(term1,term2)
        minSubsumer <- lcas
        
        if(length(lcas) >0 && !is.na(minSubsumer) && minSubsumer!="NA") {
          # get distance of minSubsumer from root
          #if (minSubsumer==rootTerm)
          #  distMinSubsumerRoot <- 0
          #else {
          #  if (length(g1@nodes[minSubsumer])>0 && length(g1@nodes[rootTerm]) >0)
          #    distMinSubsumerRoot<- sp.between(g1,minSubsumer,rootTerm,TRUE)[[1]]$length
          #  else
          #    distMinSubsumerRoot<- Inf
          #}
        
          # get distance from term1, term2 to minSubsumer
          #minSubsumer <- minSubsumerIC
          dist1 = dijkstra.sp(g3,term1)$distances[minSubsumer]
          dist2 = dijkstra.sp(g3,term2)$distances[minSubsumer]
          # get distance between the two terms
          dist = dist1 + dist2
          # get number of offsprings of minSubsumer
          numOffspringsMinSubsumer <- length(offspring[minSubsumer][[1]])  
          if (numOffspringsMinSubsumer < 1) 
            numOffspringsMinSubsumer=1
          # get IC (Resnik) of minSubsumer
          #res<-goSim(term1, term2, ont = ont, organism = "yeast", measure = "Resnik")
          #simPL = exp(-0.8/(distMinSubsumerRoot) - 0.2*dist)
          simOffspring = (1-(log(numOffspringsMinSubsumer)/log(maxOffspring)))*exp(-0.2*dist)
          #simIC = res * exp(-0.2*dist)
          } # if(!is.na(minSubsumer) && minSubsumer!="NA")
          # two terms may be disjoint
          else {
            distMinSubsumerRoot<- Inf
            simOffspring=0
          }
        # write to cache
        similarity_cache<-rbind(similarity_cache, c(term1, term2, ontology, simOffspring))
      }
      tempMatrixOffspring[i,j] = simOffspring
    } # j loop 
  } # i loop
  
  # Offspring
  #rowMax<- apply(tempMatrixOffspring, MARGIN=1, max, na.rm=TRUE)
  rowMax <- unlist(apply(tempMatrixOffspring, MARGIN=1, FUN=function(x) {if (length(which(is.na(x)))>0) x<-x[-which(is.na(x))]; if(length(x)>0) max(x) }))
  #colMax <- apply(tempMatrixOffspring, MARGIN=2, max, na.rm=TRUE)
  colMax <- unlist(apply(tempMatrixOffspring, MARGIN=2, FUN=function(x) {if (length(which(is.na(x)))>0) x<-x[-which(is.na(x))]; if(length(x)>0) max(x) }))
  overallMax <- max(as.double(rowMax), as.double(colMax))
  OffspringSimBMA = (sum(as.double(rowMax))+sum(as.double(colMax)))/(length(rowMax)+length(colMax))
  OffspringSimMax = overallMax
  
  ###### using GOSemSim package ######################
  resnik <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Resnik", combine="max")
  if (is.infinite(resnik)) resnik <-1
  resnikBMA <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Resnik", combine="BMA")
  if (is.infinite(resnikBMA)) resnikBMA <-1
  lin <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Lin", combine="max")
  if (is.infinite(lin)) lin <-1
  linBMA <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Lin", combine="BMA")
  if (is.infinite(linBMA)) linBMA <-1
  schlicker <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Rel", combine="max")
  if (is.infinite(schlicker)) schlicker <-1
  schlickerBMA <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Rel", combine="BMA")
  if (is.infinite(schlickerBMA)) schlickerBMA <-1
  jiang <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Jiang", combine="max")
  if (is.infinite(jiang)) jiang <-1
  jiangBMA <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Jiang", combine="BMA")
  if (is.infinite(jiangBMA)) jiangBMA <-1
  wang <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Wang", combine="max")
  if (is.infinite(wang)) wang <-1
  wangBMA <- mgoSim(GO1=goTerms1, GO2=goTerms2, ont=ont, organism="yeast", measure="Wang", combine="BMA")
  if (is.infinite(wangBMA)) wangBMA <-1
  ####################################################
  #intersectionGO<-intersect(goTerms1,goTerms2)
  #unionGO<-union(goTerms1,goTerms2)
  #jaccard[i,]<-c(gene1, gene2, paste(intersectionGO,collapse="/"), paste(unionGO,collapse="/"))
  

  similarity[ctr,]<-c(start+ctr, gene1, gene2,OffspringSimMax, OffspringSimBMA, resnik, resnikBMA, lin, linBMA, schlicker, schlickerBMA, jiang, jiangBMA, wang, wangBMA)
}

colnames(similarity)<-c("number","gene1","gene2", "OffspringSimMax", "OffspringSimBMA", "Resnik", "ResnikBMA", "Lin", "LinBMA", "Schlicker", "SchlickerBMA", "Jiang", "JiangBMA", "Wang", "WangBMA")
write.table(similarity,paste("similarity",ontology,"_full_start_",start,".txt",sep=""))
write.table(similarity_cache, paste("../similarity_cache_",ontology,".txt", sep=""))


