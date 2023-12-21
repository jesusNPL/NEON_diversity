library(phytools)
library(geiger)

##### Auxiliary functions #####

### Compare models using simulated traits
demonTreeFit <- function(tree, traitSims, ncores = 20) {
  
  models <- list()
  
  for(i in 1:ncol(traitSims)) { 
    
    print(paste0("Model comparison using simulated traits ", i))
    
    m_bm <- fitContinuous(phy = phy, 
                          dat = trait_sims[, i], 
                          model = "BM", 
                          ncores = ncores)
    
    m_ou <- fitContinuous(phy = phy, 
                          dat = trait_sims[, i], 
                          model = "OU", 
                          ncores = ncores)
    
    mods <- c(m_bm$opt$aic, m_ou$opt$aic) 
    names(mods) <- c("BM", "OU")
    aics <- aicw(mods)
    aics$Models <- c("BM", "OU")
    aics$Simulation <- i 
    print(aics)
    models[[i]] <- data.frame(aics)
    
  }
  
  models <- do.call(rbind, models)
  
  return(models)
  
  paste0("Model comparison finished without issues...!")
  
}

##### Load data #####
tol = 1e-9

neon_trs <- read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S2.nex")
neon_tree <- read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S3.nex")

neon_tree <- ape::multi2di(neon_tree, random = FALSE)
neon_tree <- phytools::force.ultrametric(neon_tree)
neon_tree$edge.length[neon_tree$edge.length < tol-9] <- tol

##### Simulate traits #####
trait_sims <- fastBM(neon_tree, nsim = 1000)

##### Run model comparison #####

NEON_trait_evol_sim <- demonTreeFit(tree = neon_tree, 
                                    traitSims = trait_sims, 
                                    ncores = 30)

### Save results
save(NEON_trait_evol_sim, file = "Results/Tree_accuracy/NEON_trait_evol_sim.RData")
write.csv(NEON_trait_evol_sim, file = "Results/Tree_accuracy/NEON_trait_evol_sim.csv")
