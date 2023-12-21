
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

tips <- data.frame(sciname = tr$tip.label) %>% 
  mutate(Taxa = sciname) %>% 
  separate(sciname, into = c("genus", "epithet"), sep = "_") %>% 
  select(Taxa, genus, epithet)
  
##### Match trait - phylogeny #####
traits <- X[X$Taxa %in% tr$tip.label, ]

traits <- full_join(traits, tips, by = c("Taxa", "genus", "epithet"))

rownames(traits) <- traits$Taxa

tr <- drop.tip(tr, setdiff(tr$tip.label, traits$Taxa))


traits2 <- traits[, 4:35]

##### Run imputation ######

pvr <- PVR::PVRdecomp(tr)

# Extract the PVRs
pvrs <- pvr@Eigen$vectors

# Combine traits and PVRs
traits_pvrs <- cbind(traits2, pvrs)


###### Imputation using missForest #####
# (note that this function have other arguments, see details in Stekhoven and Buehlmann 2012)

library(doParallel)

registerDoParallel(cores = 30)

imp <- missForest::missForest(xmis = traits_pvrs, 
                              maxiter = 100, 
                              ntree = 250, 
                              variablewise = TRUE, 
                              parallelize = "forest", 
                              verbose = TRUE)

save(imp, file = "DATA/Traits/new_Imputation_FULL.RData") 

##### Post-processing #####
load("DATA/Traits/new_Imputation_FULL.RData")

# Imputed traits
imputedTraits <- imp$ximp[, 1:32] 

# OOBerror
errors <- imp$OOBerror[1:32]
names(errors) <- names(imputedTraits)

errors 

## Extract imputed traits 
# Inits
nums1 <- seq(1, 31, by = 2)
nums2 <- seq(2, 32, by = 2)

# Mean trait imputed values
trt_mean_names <- names(imputedTraits)[nums1]
# SD trait imputed values
trt_sd_names <- names(imputedTraits)[nums2]

### Select traits 
impTraits <- imputedTraits %>% 
  select(trt_mean_names) %>% 
  mutate(species = rownames(imputedTraits)) %>% 
  select(species, everything())

impTraitsSD <- imputedTraits %>% 
  select(trt_sd_names) %>% 
  mutate(species = rownames(imputedTraits)) %>% 
  select(species, everything())

### Save trait imputation
write.csv(impTraits, file = "DATA/Traits/GapFilling/imputedTRAITS.csv")

write.csv(impTraitsSD, file = "DATA/Traits/GapFilling/imputedTRAITS_SD.csv")


