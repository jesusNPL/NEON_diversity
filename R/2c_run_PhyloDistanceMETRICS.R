library(tidyverse)
library(phytools)

### Auxiliary function
source("R/NEON_diversity/R/Functions/AlphaSpecDIV.R")

### Load data
load("DATA/matchPhyloComm/matchPhyloComm.RData")

### Inits
# Site names
sites <- names(matched_PhyComm) 
# Number of sites
nSites <- length(sites)
# Diversity order
QS <- c(0, 1, 2) # different q values

# Store results
phylo_NEON_dis_q0 <- list()
phylo_NEON_dis_q1 <- list()
phylo_NEON_dis_q2 <- list()

### Start calculations

for(j in 1:nSites) { 
  
  # Select site j
  matched <- matched_PhyComm[[j]]
  # Get plot names per site
  comNames <- rownames(matched$comm)
  
  ### Q 0
  print(paste0("Starting calculations for site ", sites[j], " under Q0"))
  
  phylo_dis_q0 <- demon_PhyloDistance(comm = matched$comm, 
                                      phylo = matched$phy, 
                                      abundance = TRUE, 
                                      Q = QS[1], 
                                      siteName = sites[j], 
                                      plotNames = comNames) 
  
  phylo_NEON_dis_q0[[j]] <- phylo_dis_q0
  
  ### Q 1
  print(paste0("Starting calculations for site ", sites[j], " under Q1"))
  
  phylo_dis_q1 <- demon_PhyloDistance(comm = matched$comm, 
                                      phylo = matched$phy, 
                                      abundance = TRUE, 
                                      Q = QS[2], 
                                      siteName = sites[j], 
                                      plotNames = comNames)  
  
  phylo_NEON_dis_q1[[j]] <- phylo_dis_q1
  
  ### Q 2
  print(paste0("Starting calculations for site ", sites[j], " under Q2"))
  
  phylo_dis_q2 <- demon_PhyloDistance(comm = matched$comm, 
                                      phylo = matched$phy, 
                                      abundance = TRUE, 
                                      Q = QS[3], 
                                      siteName = sites[j], 
                                      plotNames = comNames)  
  
  phylo_NEON_dis_q2[[j]] <- phylo_dis_q2

}

### Combine results per site
phylo_NEON_table_q0 <- do.call(rbind, phylo_NEON_dis_q0)
phylo_NEON_table_q1 <- do.call(rbind, phylo_NEON_dis_q1)
phylo_NEON_table_q2 <- do.call(rbind, phylo_NEON_dis_q2)

### Save results
save(phylo_NEON_table_q0, 
     phylo_NEON_table_q1, 
     phylo_NEON_table_q2, 
     file = "Results/phyloDiversity/Reanalyses/phylo_DIS_tables.RData")

save(phylo_NEON_dis_q0, 
     phylo_NEON_dis_q1, 
     phylo_NEON_dis_q2, 
     file = "Results/phyloDiversity/Reanalyses/phylo_DIS_lists.RData")
