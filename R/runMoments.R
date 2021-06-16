library(dplyr)
library(tidyr)

source("R/NEON_diversity/R/Functions/demon_moments.R")

#### Prepare data #####
load("DATA/RData/HARV_SS.RData") # Spectral species

### Clean vegetation community data 
load("DATA/Harvard_NEON/Data/harv_matchData.RData") # HARV data
rm(harv_taxonomy, matchTrait, matchPhylo)

harv_comm <- harv_comm[harv_comm$plotID_SPEC %in% unique(ss_threshold$plotID), ]
rownames(harv_comm) <- harv_comm$plotID_SPEC

harv_comm_clean <- harv_comm_clean[harv_comm_clean$plotID_SPEC %in% unique(ss_threshold$plotID), ]
rownames(harv_comm_clean) <- harv_comm_clean$plotID_SPEC 

##### Spectral moments #####
## Spectra data range from column 4 to 430
plotNames <- unique(ss_threshold$plotID)
nPlots <- length(plotNames)
Range <- 5:430 ## Spectra data range from column 5 to 430

moment_spec <- demon_momentsSPEC(spectra = ss_threshold, nPlots = nPlots, 
                       plotNames = plotNames, specRange = Range)

##### Dimension moments #####
harv_comm_clean[1:5, 1:5]

phy <- ape::read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S3.nex")
phyHarv <- ape::drop.tip(phy, setdiff(phy$tip.label, names(harv_comm_clean)))

traits <- read.csv("DATA/Traits/Traits_BIEN/GapFilling/imputedTRAITS.csv")
traitHarv <- traits[traits$Species %in% phyHarv$tip.label, ]

trtHarv <- data.frame(traits[, c("Species", "mean_LAreaPLDryMass")])

moment_div <- demon_momentsTraitPhylo(comm = harv_comm_clean[, 4:121], 
                                      trait = trtHarv, 
                                      phylo = phyHarv, 
                                      plotNames = plotNames, 
                                      nPlots = nPlots)

### combine results
moments_tbl <- makeTable_moments(divTaxo = moment_div$moments_taxo, 
                                 divTrait = moment_div$moments_traits, 
                                 divPhylo = moment_div$moments_phylo, 
                                 divSpec = moment_spec)

##### Make regressions #####
demon_BayCorrMOMENTS(momentsRES = moments_tbl, 
                     nChains = 4, 
                     nIters = 2000, 
                     nCores = 10, 
                     pathSaveTaxo = "Results/HARV/Bay_taxoSpec_moments_correlations.RData", 
                     pathSaveTrait = "Results/HARV/Bay_traitSpec_moments_correlations.RData", 
                     pathSavePhylo = "Results/HARV/Bay_phyloSpec_moments_correlations.RData")

