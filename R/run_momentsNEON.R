library(tidyverse)

source("R/NEON_diversity/R/Functions/demon_moments.R")

### Load data
load("DATA/matchPhyloComm/matchPhyloComm.RData")
load("DATA/matchTraitComm/matchTraitComm.RData")

### Phylogeny
phyNEON <- ape::read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S3.nex")

### Inits
sites <- names(matched_PhyComm)
nSites <- length(sites)

moments_NEON <- list()

for(i in 1:nSites) { 
  
  print(paste0("Run moments for site ", sites[i], " ..."))
  
  comTMP <- matchedTrait[[2]][[i]]
  traitTMP <- matchedTrait[[1]][[i]]
  
  plotNames <- rownames(comTMP)
  nPlots <- length(plotNames)
  
  momSITE <- demon_momentsTraitPhylo(comm = comTMP, 
                                     trait = traitTMP, 
                                     phylo = phyNEON, 
                                     plotNames = plotNames, 
                                     nPlots = nPlots, 
                                     siteName = sites[i], 
                                     maxIter = 100)

  #momSITE$Site <- sites[i]
  moments_NEON[[i]] <- momSITE
  
}


save(moments_NEON, file = "Results/moments/moments_NEON.RData")
