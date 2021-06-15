
## Load packages
require(geiger)
require(PVR)
require(phytools)
require(missForest)

source("R/NEON_diversity/R/Functions/demonImputation.R")

X <- readRDS("DATA/Traits/Traits_BIEN/traits_BIEN_NEON_spp_genus_combined.rds")
#Y <- readRDS("DATA/Traits/Traits_BIEN/spp_traits_BIEN_TPL.rds")
tr <- ape::read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S3.nex")
trs <- read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S2.nex")

##### Match trait - phylogeny #####
traits <- X[X$Taxa %in% tr$tip.label, ]
rownames(traits) <- traits$Taxa

tr <- drop.tip(tr, setdiff(tr$tip.label, traits$Taxa))
trs <- lapply(trs, drop.tip, setdiff(trs[[1]]$tip.label, traits$Taxa))
class(trs) <- "multiPhylo"
#
traits2 <- traits[, 4:35]
#traitSel <- traits[, 4:5]

##### Run imputation ######

pvrRES <- makePVR(phylo = tr, multiPhylo = trs)

# Inits
nums1 <- seq(1, 31, by = 2)
nums2 <- seq(2, 32, by = 2)

traitImputation <- list()

for(j in 1:16) {
  
  sel <- names(traits2)[nums1[j]:nums2[j]]
  
  print(paste0("Run trait imputation for ", sel[1], " ..."))
  #print(sel[1])
  
  traitSel <- traits2[, sel]
  
  traitImputation[[j]] <- demon_ImpTrait(traits = traitSel, traitNames = names(traitSel), 
                                         phylo = tr, multiPhylo = trs, 
                                         iters = 100, ntrees = 250, 
                                         pvr = pvrRES[[1]], pvrSamp = pvrRES[[2]])
  
}

## Save results
save(traitImputation, file = "DATA/Traits/Traits_BIEN/GapFilling/Imputation_Traits.RData")

##### Post-processing imputation #####

load("DATA/Traits/Traits_BIEN/GapFilling/Imputation_Traits.RData")

## Extract imputed traits from best tree
nTraits = 16

impTraits <- list()
impTraitSD <- list()

for(i in 1:nTraits) {
  impTraits[[i]] <- traitImputation[[i]]$impTraits[1]
  impTraitSD[[i]] <- traitImputation[[i]]$impTraits[2]
  
}

impTraits <- do.call(cbind, impTraits)
write.csv(impTraits, file = "DATA/Traits/Traits_BIEN/GapFilling/imputedTRAITS.csv")

impTraitSD <- do.call(cbind, impTraitSD)
write.csv(impTraitSD, file = "DATA/Traits/Traits_BIEN/GapFilling/imputedTRAITS_SD.csv")

## Extract imputed traits from sample of trees 
nums1 <- seq(1, 200, by = 2)
nums2 <- seq(2, 200, by = 2)

impTraits_100 <- list()
impTraitSD_100 <- list()

for(i in 1:nTraits) {
  tmp <- traitImputation[[i]]$impMulti 
  
  # Mean values
  X <- do.call(cbind, tmp)[nums1]
  X <- data.frame(rowMeans(X))
  names(X) <- names(tmp[[1]])[1] 
  impTraits_100[[i]] <- X 
  # SD values
  Y <- do.call(cbind, tmp)[nums2]
  Y <- data.frame(rowMeans(Y))
  names(Y) <- names(tmp[[1]])[2] 
  impTraitSD_100[[i]] <- Y 

}

impTraits_100 <- do.call(cbind, impTraits_100)
impTraitSD_100 <- do.call(cbind, impTraitSD_100)




