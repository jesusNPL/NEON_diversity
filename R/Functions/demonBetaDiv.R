##### Function to run multiple dimensions of beta diversity ##### 
demonBetaDiv <- function(PhyComData, TraitComData) { 
  
  ### Taxonomic beta diversity 
  # Abundance
  taxBetaAbun <- list()
  taxBetaAbun_bal <- list()
  taxBetaAbun_gra <- list()
  # Presence-Absence
  taxBetaPA <- list()
  taxBetaPA_sim <- list()
  taxBetaPA_sne <- list()
  
  ### Phylogenetic beta diversity
  phyBeta <- list()
  phyBeta_sim <- list()
  phyBeta_sne <- list()
  
  ### Trait beta diversity
  traitBeta <- list()
  traitBeta_sim <- list()
  traitBeta_sne <- list()
  
  ### Inits
  siteName <- names(PhyComData)
  nSites <- length(siteName)
  
  ### Run
  for(i in 1:nSites) { 
    
    print(paste0("Taxonomic beta diversity for site ", siteName[i]))
    
    ##### Taxonomic beta diversity #####
    comAB <- PhyComData[[i]]$comm
    comPA <- comAB
    comPA[comPA > 0] <- 1 
    
    ## Abundance beta diversity
    betaAB <- betapart::beta.pair.abund(comAB, index.family = "bray")
    
    bray <- betaAB$beta.bray
    bray <- as.matrix(bray)
    bray <- data.frame(metricVal = apply(bray, 2, FUN = mean))
    bray$type <- "bray_curtis"
    bray$plotID <- rownames(bray)
    bray$Site <- siteName[i]
    
    bal <- betaAB$beta.bray.bal 
    bal <- as.matrix(bal)
    bal <- data.frame(metricVal = apply(bal, 2, FUN = mean))
    bal$type <- "bray_bal"
    bal$plotID <- rownames(bal)
    bal$Site <- siteName[i]
    
    gra <- betaAB$beta.bray.gra
    gra <- as.matrix(gra)
    gra <- data.frame(metricVal = apply(gra, 2, FUN = mean))
    gra$type <- "bray_gra"
    gra$plotID <- rownames(gra)
    gra$Site <- siteName[i]
    
    # Store results
    taxBetaAbun[[i]] <- bray
    taxBetaAbun_bal[[i]] <- bal
    taxBetaAbun_gra[[i]] <- gra
    
    ### PA beta diversity
    betaPA <- betapart::beta.pair(comPA, index.family = "sorensen")
    
    sor <- betaPA$beta.sor
    sor <- as.matrix(sor)
    sor <- data.frame(metricVal = apply(sor, 2, FUN = mean))
    sor$type <- "sorensen"
    sor$plotID <- rownames(sor)
    sor$Site <- siteName[i]
    
    sim <- betaPA$beta.sim 
    sim <- as.matrix(sim)
    sim <- data.frame(metricVal = apply(sim, 2, FUN = mean))
    sim$type <- "turnover"
    sim$plotID <- rownames(sim)
    sim$Site <- siteName[i]
    
    sne <- betaPA$beta.sne
    sne <- as.matrix(sne)
    sne <- data.frame(metricVal = apply(sne, 2, FUN = mean))
    sne$type <- "nestedness" 
    sne$plotID <- rownames(sne)
    sne$Site <- siteName[i]
    
    # Store results
    taxBetaPA[[i]] <- sor
    taxBetaPA_sim[[i]] <- sim
    taxBetaPA_sne[[i]] <- sne 
    
    ##### Phylogenetic beta diversity ##### 
    print(paste0("Phylogenetic beta diversity for site ", siteName[i]))
    
    comAB <- PhyComData[[i]]$comm
    comPA <- comAB
    comPA[comPA > 0] <- 1 
    
    phy <- PhyComData[[i]]$phy
    
    phybetaPA <- betapart::phylo.beta.pair(x = comPA, tree = phy, index.family = "sorensen")
    
    physor <- phybetaPA$phylo.beta.sor
    physor <- as.matrix(physor) 
    physor <- data.frame(metricVal = apply(physor, 2, FUN = mean))
    physor$type <- "sorensen"
    physor$plotID <- rownames(physor)
    physor$Site <- siteName[i]
    
    physim <- phybetaPA$phylo.beta.sim 
    physim <- as.matrix(physim)
    physim <- data.frame(metricVal = apply(physim, 2, FUN = mean))
    physim$type <- "turnover"
    physim$plotID <- rownames(physim)
    physim$Site <- siteName[i]
    
    physne <- phybetaPA$phylo.beta.sne
    physne <- as.matrix(physne) 
    physne <- data.frame(metricVal = apply(physne, 2, FUN = mean))
    physne$type <- "nestedness"
    physne$plotID <- rownames(physne)
    physne$Site <- siteName[i]
    
    # Store results 
    phyBeta[[i]] <- physor
    phyBeta_sim[[i]] <- physim
    phyBeta_sne[[i]] <- physne
    
    ##### Trait beta diversity #####
    print(paste0("Trait beta diversity for site ", siteName[i]))
    
    comAB <- TraitComData[[i]]$comm
    comPA <- comAB
    comPA[comPA > 0] <- 1 
    
    phy <- TraitComData[[i]]$phy
    
    trtbetaPA <- betapart::phylo.beta.pair(x = comPA, tree = phy, index.family = "sorensen")
    
    trtsor <- trtbetaPA$phylo.beta.sor
    trtsor <- as.matrix(trtsor)
    trtsor <- data.frame(metricVal = apply(trtsor, 2, FUN = mean))
    trtsor$type <- "sorensen"
    trtsor$plotID <- rownames(trtsor)
    trtsor$Site <- siteName[i]
    
    trtsim <- trtbetaPA$phylo.beta.sim 
    trtsim <- as.matrix(trtsim) 
    trtsim <- data.frame(metricVal = apply(trtsim, 2, FUN = mean))
    trtsim$type <- "turnover"
    trtsim$plotID <- rownames(trtsim)
    trtsim$Site <- siteName[i]
    
    trtsne <- trtbetaPA$phylo.beta.sne
    trtsne <- as.matrix(trtsne) 
    trtsne <- data.frame(metricVal = apply(trtsne, 2, FUN = mean))
    trtsne$type <- "nestedness"
    trtsne$plotID <- rownames(trtsne)
    trtsne$Site <- siteName[i]
    
    # Store results 
    traitBeta[[i]] <- trtsor
    traitBeta_sim[[i]] <- trtsim
    traitBeta_sne[[i]] <- trtsne
    
  } 
  
  ##### Merge results #####
  ## Taxonomic abundance
  taxBetaAbun <- do.call(rbind, taxBetaAbun)
  taxBetaAbun_bal <- do.call(rbind, taxBetaAbun_bal)
  taxBetaAbun_gra <- do.call(rbind, taxBetaAbun_gra)
  
  taxBetaAbundance <- rbind(taxBetaAbun, taxBetaAbun_bal, taxBetaAbun_gra)
  rownames(taxBetaAbundance) <- NULL
  
  ## Taxonomic presence-absence
  taxBetaPA <- do.call(rbind, taxBetaPA)
  taxBetaPA_sim <- do.call(rbind, taxBetaPA_sim)
  taxBetaPA_sne <- do.call(rbind, taxBetaPA_sne)
  
  taxBetaPresAbs <- rbind(taxBetaPA, taxBetaPA_sim, taxBetaPA_sne)
  rownames(taxBetaPresAbs) <- NULL
  
  ## Phylogenetic
  phyBeta <- do.call(rbind, phyBeta)
  phyBeta_sim <- do.call(rbind, phyBeta_sim)
  phyBeta_sne <- do.call(rbind, phyBeta_sne)
  
  phyloBetaDiv <- rbind(phyBeta, phyBeta_sim, phyBeta_sne)
  rownames(phyloBetaDiv) <- NULL
  
  ## Trait
  traitBeta <- do.call(rbind, traitBeta)
  traitBeta_sim <- do.call(rbind, traitBeta_sim) 
  traitBeta_sne <- do.call(rbind, traitBeta_sne)
  
  traitBetaDiv <- rbind(traitBeta, traitBeta_sim, traitBeta_sne) 
  rownames(traitBetaDiv) <- NULL
  
  ##### return results #####
  res <- list("taxonomicAbundance" = taxBetaAbundance, 
              "taxonomicPresAbs" = taxBetaPresAbs, 
              "phylogenetic" = phyloBetaDiv, 
              "trait" = traitBetaDiv)
  
  return(res) 
  
}

##### Load and run beta diversity #####
#load("output/BetaDiv/data_betaDiv.RData")

#betaDiv_NEON <- demonBetaDiv(PhyComData = ComPhy_data, 
 #                            TraitComData = ComTrt_data)

#save(betaDiv_NEON, file = "output/BetaDiv/results_betaDiv.RData")

##### Function to run multiple dimensions of beta diversity multisite ##### 
demonBetaDivMulti <- function(PhyComData, TraitComData) { 
  
  metrics <- data.frame(matrix(ncol = 15, nrow = length(PhyComData)))
  names(metrics) <- c("Site", "meanSR", "meanRangesize",
                      "bray_curtis", "bray_curtis_bal", "bray_curtis_gra", 
                      "taxoSorensen", "taxoTurnover", "taxoNestedness", 
                      "phyloSorensen", "phyloTurnover", "phyloNestedness", 
                      "traitSorensen", "traitTurnover", "traitNestedness")
  
  ### Inits
  siteName <- names(PhyComData)
  nSites <- length(siteName)
  
  ### Run
  for(i in 1:nSites) { 
    
    metrics[i, 1] <- siteName[i]
    
    print(paste0("Taxonomic beta diversity for site ", siteName[i]))
    
    ##### Taxonomic beta diversity #####
    comAB <- PhyComData[[i]]$comm
    comPA <- comAB
    comPA[comPA > 0] <- 1 
    
    metrics[i, 2] <- mean(rowSums(comPA))
    metrics[i, 3] <- mean(colSums(comPA))
    ## Abundance beta diversity
    betaAB <- betapart::beta.multi.abund(comAB, index.family = "bray")
    
    metrics[i, 4] <- betaAB$beta.BRAY
    metrics[i, 5] <- betaAB$beta.BRAY.BAL
    metrics[i, 6] <- betaAB$beta.BRAY.GRA
    
    ### PA beta diversity
    betaPA <- betapart::beta.multi(comPA, index.family = "sorensen")
    
    metrics[i, 7] <- betaPA$beta.SOR
    metrics[i, 8] <- betaPA$beta.SIM
    metrics[i, 9] <- betaPA$beta.SNE
    
    ##### Phylogenetic beta diversity ##### 
    print(paste0("Phylogenetic beta diversity for site ", siteName[i]))
    
    comAB <- PhyComData[[i]]$comm
    comPA <- comAB
    comPA[comPA > 0] <- 1 
    
    phy <- PhyComData[[i]]$phy
    
    phybetaPA <- betapart::phylo.beta.multi(x = comPA, tree = phy, index.family = "sorensen")
    
    metrics[i, 10] <- phybetaPA$phylo.beta.SOR
    metrics[i, 11] <- phybetaPA$phylo.beta.SIM
    metrics[i, 12] <- phybetaPA$phylo.beta.SNE
    
    ##### Trait beta diversity #####
    print(paste0("Trait beta diversity for site ", siteName[i]))
    
    comAB <- TraitComData[[i]]$comm
    comPA <- comAB
    comPA[comPA > 0] <- 1 
    
    phy <- TraitComData[[i]]$phy
    
    trtbetaPA <- betapart::phylo.beta.multi(x = comPA, tree = phy, index.family = "sorensen")
    
    metrics[i, 13] <- trtbetaPA$phylo.beta.SOR
    metrics[i, 14] <- trtbetaPA$phylo.beta.SIM
    metrics[i, 15] <- trtbetaPA$phylo.beta.SNE
    
  } 
  
  return(metrics) 
  
}

##### Load and run beta diversity #####
#load("output/BetaDiv/data_betaDiv.RData")

#betaDivMulti_NEON <- demonBetaDivMulti(PhyComData = ComPhy_data, 
 #                                      TraitComData = ComTrt_data)

#write.csv(betaDivMulti_NEON, file = "output/BetaDiv/results_betaDivMulti.csv")
