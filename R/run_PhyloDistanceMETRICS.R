library(tidyverse)
library(picante)
library(fundiversity)
library(phytools)
source("R/NEON_diversity/R/Functions/AlphaSpecDIV.R")

### Load data
load("DATA/matchPhyloComm/matchPhyloComm.RData")

### Inits
sites <- names(matched_PhyComm)
nSites <- length(sites)
QS <- c(0, 1, 2, 3) # different q values

phylo_NEON_dis_q0 <- list()
phylo_NEON_dis_q1 <- list()
phylo_NEON_dis_q2 <- list()
phylo_NEON_dis_q3 <- list()

for(j in 1:nSites) { 
  
  matched <- matched_PhyComm[[j]]
  
  comNames <- rownames(matched$comm)
  
  print(paste0("Starting calculations for site ", sites[j], "..."))
  
  ### Q 0
  phylo_dis_q0 <- demon_PhyloDistance(comm = matched$comm, 
                                             phylo = matched$phy, 
                                             abundance = TRUE, Q = QS[1], 
                                             plotNames = comNames)  
  phylo_dis_q0$Site <- sites[j]
  phylo_NEON_dis_q0[[j]] <- phylo_dis_q0
  
  ### Q 1
  phylo_dis_q1 <- demon_PhyloDistance(comm = matched$comm, 
                                      phylo = matched$phy, 
                                      abundance = TRUE, Q = QS[2], 
                                      plotNames = comNames)  
  phylo_dis_q1$Site <- sites[j]
  phylo_NEON_dis_q1[[j]] <- phylo_dis_q1
  
  ### Q 2
  phylo_dis_q2 <- demon_PhyloDistance(comm = matched$comm, 
                                      phylo = matched$phy, 
                                      abundance = TRUE, Q = QS[3], 
                                      plotNames = comNames)  
  phylo_dis_q2$Site <- sites[[j]]
  phylo_NEON_dis_q2[[j]] <- phylo_dis_q2
  
  ###Q 3
  phylo_dis_q3 <- demon_PhyloDistance(comm = matched$comm, 
                                      phylo = matched$phy, 
                                      abundance = TRUE, Q = QS[4], 
                                      plotNames = comNames)  
  phylo_dis_q3$Site <- sites[j]
  phylo_NEON_dis_q3[[j]] <- phylo_dis_q3
  

}


for(i in 1:nSites) {
  phylo_NEON_dis_q0[[i]]$Site <- sites[i]
  phylo_NEON_dis_q1[[i]]$Site <- sites[i]
  phylo_NEON_dis_q2[[i]]$Site <- sites[i]
  phylo_NEON_dis_q3[[i]]$Site <- sites[i]
  
}

phylo_NEON_table_q0 <- do.call(rbind, phylo_NEON_dis_q0)
phylo_NEON_table_q1 <- do.call(rbind, phylo_NEON_dis_q1)
phylo_NEON_table_q2 <- do.call(rbind, phylo_NEON_dis_q2)
phylo_NEON_table_q3 <- do.call(rbind, phylo_NEON_dis_q3)

save(phylo_NEON_table_q0, phylo_NEON_table_q1, 
     phylo_NEON_table_q2, phylo_NEON_table_q3, 
     file = "Results/phyloDiversity/phylo_DIS_tables.RData")

save(phylo_NEON_dis_q0, phylo_NEON_dis_q1, 
     phylo_NEON_dis_q2, phylo_NEON_dis_q3, 
     file = "Results/phyloDiversity/phylo_DIS_lists.RData")
