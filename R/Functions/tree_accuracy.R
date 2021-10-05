### Compare models using simulated traits
demonTreeAdequacy <- function(tree, nSims, ncores = 20, brlens_tol = 1e-9) {
  
  ### Correcting tree
  phy <- ape::multi2di(tree, random = FALSE)
  phy <- phytools::force.ultrametric(phy)
  phy$edge.length[phy$edge.length < brlens_tol-9] <- brlens_tol
  ape::ladderize(phy)
  
  ### Simulate trait values
  trait_sims <- phytools::fastBM(phy, nsim = nSims)
  # remove NAs caused for some very short branch lenghts
  trait_sims <- na.omit(trait_sims)
  
  spp_keep <- rownames(trait_sims)
  # Prune phylogeny
  phy <- ape::drop.tip(phy, setdiff(phy$tip.label, spp_keep))
  
  models <- list()
  
  for(i in 1:nSims) { 
    print(paste0("Model comparison using simulated trait values ", i))
    
    m_bm <- geiger::fitContinuous(phy = phy, 
                          dat = trait_sims[, i], 
                          model = "BM", 
                          ncores = ncores)
    
    m_ou <- geiger::fitContinuous(phy = phy, 
                          dat = trait_sims[, i], 
                          model = "OU", 
                          ncores = ncores)
    
    mods <- c(m_bm$opt$aic, m_ou$opt$aic) 
    names(mods) <- c("BM", "OU")
    aics <- geiger::aicw(mods)
    aics$Models <- c("BM", "OU")
    aics$Simulation <- i 
    print(aics)
    models[[i]] <- aics
    
  }
  
  models <- do.call(rbind, models)
  return(models)
  print(paste0("Model comparison finished without issues...!"))
}

#test <- demonTreeAdequacy(tree = neon_tree, nSims = 5, ncores = 4)
