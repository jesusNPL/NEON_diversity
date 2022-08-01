
## Load packages
require(geiger)
require(PVR)
require(phytools)
require(missForest)
library(Taxonstand)

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

### Extract errors
nums1 <- seq(1, 31, by = 2)

errors <- numeric(length = nTraits)

for(i in 1:nTraits) { 
  errors[i] <- traitImputation[[i]]$NRMSE
}

names(errors) <- names(traits2[, c(nums1)])

errors

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


### Extract errors
nums1 <- seq(1, 31, by = 2)

errors_Samp <- list()

for(i in 1:nTraits) { 
  errors_Samp[[i]] <- traitImputation[[i]]$NRMSEs
}

names(errors) <- names(traits2[, c(nums1)])

errors

##### Re-run imputation using mean values #####
library(tidyverse)
library(missForest)
library(geiger)

X <- readRDS("DATA/Traits/Traits_BIEN/traits_BIEN_NEON_spp_genus_combined.rds")
#Y <- readRDS("DATA/Traits/Traits_BIEN/spp_traits_BIEN_TPL.rds")
tr <- ape::read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S3.nex")

trt_world <- read.csv("../Lucie_Jesus/RCN_evolutionary_legacy/Data/TryData/FinalDATA/traitsMATCH_checkNOnas_FINAL.csv")

##### Match trait - phylogeny #####
traits <- X[X$Taxa %in% tr$tip.label, ] 

trait_sel <- trt_world[trt_world$Species %in% unique(traits$Taxa), ]

lpm <- trt_world %>% 
  select(Species, lpm) %>% 
  drop_na() %>% 
  group_by(Species) %>% 
  summarize(mean_lmp = mean(lpm))

sla <- trt_world %>% 
  select(Species, sla) %>% 
  drop_na() %>% 
  group_by(Species) %>% 
  summarize(mean_sla = mean(sla))

lnm <- trt_world %>% 
  select(Species, lnm) %>% 
  drop_na() %>% 
  group_by(Species) %>% 
  summarize(mean_lnm = mean(lnm))

lpm <- data.frame(lpm)
sla <- data.frame(sla)
lnm <- data.frame(lnm)

traits <- full_join(X, lpm, by = c("Taxa" = "Species"))
traits <- full_join(traits, sla, by = c("Taxa" = "Species"))
traits <- full_join(traits, lnm, by = c("Taxa" = "Species"))

traits <- full_join(traits, lpm, by = c("Taxa" = "Species"))
rownames(traits) <- traits$Taxa

traits <- traits[traits$Taxa %in% tr$tip.label, ]


tr <- drop.tip(tr, setdiff(tr$tip.label, traits$Taxa))
#
traits2 <- traits[, 4:36]
nums1 <- seq(1, 33, by = 2)

traits2 <- traits2[, nums1] 
traits2 <- apply(traits2, 2, FUN = scale) 
rownames(traits2) <- traits$Taxa
#traitSel <- traits[, 4:5]

##### Run imputation ######
## Imputation
# Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
pvrRES <- PVR::PVRdecomp(tr)

# Extract the PVRs
pvrs <- pvrRES@Eigen$vectors

# Combine traits and PVRs
traits_pvrs <- cbind(traits2, pvrs)

# Imputation using missForest (note that this function have other arguments, see details in Stekhoven and Buehlmann 2012)
library(doParallel)
registerDoParallel(cores = 28)

imp <- missForest::missForest(xmis = traits_pvrs, 
                              maxiter = 100, 
                              ntree = 250, 
                              variablewise = TRUE, 
                              parallelize = "forest", 
                              verbose = FALSE)

save(imp, file = "DATA/Traits/new_Imputation.RData") 


impNO <- missForest::missForest(xmis = traits2, 
                              maxiter = 100, 
                              ntree = 250, 
                              variablewise = TRUE, 
                              parallelize = "forest", 
                              verbose = FALSE)

save(impNO, file = "DATA/Traits/new_Imputation_NOphy.RData")


impTraits <- imp$ximp[, traitNames]
# OOBerror
err <- imp$OOBerror

imp2 <- missForest::missForest(xmis = traits_pvrs, 
                              maxiter = 100, 
                              ntree = 250, 
                              variablewise = FALSE, 
                              parallelize = "forest", 
                              verbose = FALSE)

save(imp2, file = "DATA/Traits/new_Imputation2.RData")
