##### Run functional diversity #####
library(tidyverse)
library(picante)
library(fundiversity)
library(phytools)
library(gawdis)

source("R/NEON_diversity/R/Functions/AlphaSpecDIV.R")
source("R/NEON_diversity/R/Functions/auxiliar.R")

### load data
load("DATA/matchTraitComm/matchTraitComm.RData")

### inits 
sites <- names(matchedTrait$traitNEON)
nSites <- length(sites)
QS <- c(0, 1, 2, 3) # different q values


trait_NEON_dis_q0 <- list()
trait_NEON_dis_q1 <- list()
trait_NEON_dis_q2 <- list()
trait_NEON_dis_q3 <- list()

for(i in 1:nSites) {
  
  ### Inits
  comTMP <- matchedTrait[[2]][[i]]
  traitTPM <- matchedTrait[[1]][[i]]
  rownames(traitTPM) <- traitTPM$Species
  
  plotNames <- rownames(comTMP)
  nPlots <- length(plotNames)
  
  print(paste0("Starting calculations under for site ", sites[i], " under Q0"))
  
  ## Q 0
  trait_NEON_dis_q0[[i]] <- demon_TraitDistance(comm = comTMP, 
                                                trait = traitTPM[, 2:ncol(traitTPM)], 
                                                abundance = TRUE, 
                                                Q = QS[1], plotNames = plotNames, 
                                                nPlots = nPlots, gowDist = TRUE, gawDist = FALSE) 
  
  print(paste0("Starting calculations under for site ", sites[i], " under Q1"))
  
  ## Q 1
  trait_NEON_dis_q1[[i]] <- demon_TraitDistance(comm = comTMP, 
                                                trait = traitTPM[, 2:ncol(traitTPM)], 
                                                abundance = TRUE, 
                                                Q = QS[2], plotNames = plotNames, 
                                                nPlots = nPlots, gowDist = TRUE, gawDist = FALSE) 
  
  print(paste0("Starting calculations under for site ", sites[i], " under Q2"))
  
  ## Q 2
  trait_NEON_dis_q2[[i]] <- demon_TraitDistance(comm = comTMP, 
                                                trait = traitTPM[, 2:ncol(traitTPM)], 
                                                abundance = TRUE, 
                                                Q = QS[3], plotNames = plotNames, 
                                                nPlots = nPlots, gowDist = TRUE, gawDist = FALSE) 
  
  print(paste0("Starting calculations under for site ", sites[i], " under Q3"))
  
  ## Q 3
  trait_NEON_dis_q3[[i]] <- demon_TraitDistance(comm = comTMP, 
                                                trait = traitTPM[, 2:ncol(traitTPM)], 
                                                abundance = TRUE, 
                                                Q = QS[4], plotNames = plotNames, 
                                                nPlots = nPlots, gowDist = TRUE, gawDist = FALSE)
}

### Add sites to each calculation 

for(i in 1:nSites) {
  trait_NEON_dis_q0[[i]]$Site <- sites[i]
  trait_NEON_dis_q1[[i]]$Site <- sites[i]
  trait_NEON_dis_q2[[i]]$Site <- sites[i]
  trait_NEON_dis_q3[[i]]$Site <- sites[i]
}

trait_NEON_dis_table_q0 <- do.call(rbind, trait_NEON_dis_q0)
trait_NEON_dis_table_q1 <- do.call(rbind, trait_NEON_dis_q1)
trait_NEON_dis_table_q2 <- do.call(rbind, trait_NEON_dis_q2)
trait_NEON_dis_table_q3 <- do.call(rbind, trait_NEON_dis_q3)

### Save results 
save(trait_NEON_dis_table_q0, trait_NEON_dis_table_q1, 
     trait_NEON_dis_table_q2, trait_NEON_dis_table_q3, 
     file = "Results/traitDiversity/trait_DIS_tables.RData")

save(trait_NEON_dis_q0, trait_NEON_dis_q1, 
     trait_NEON_dis_q2, trait_NEON_dis_q3, 
     file = "Results/traitDiversity/trait_DIS_lists.RData")

