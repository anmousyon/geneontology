while(length(goTerms)>0){
  # always start with first term and then remove when done
  #t<-goTerms[8]
  t<-"GO:0006338"
  if (length(goTerms)>1) {
    otherTerms <- goTerms[-which(goTerms==t)]
    lca<-as.character(sapply(otherTerms, FUN=function(k) findLCA(t,k)))
    lcadepth<-sapply(lca, findDepth)
    lcaOffspringCount<-sapply(lca, findOffspringCount)
    # find lca with max depth
    if (max(lcadepth) > depthThreshold) {
      browser()
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