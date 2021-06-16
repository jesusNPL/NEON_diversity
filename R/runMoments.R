library(dplyr)
library(tidyr)
source("R/NEON_diversity/R/Functions/demon_moments.R")

load("DATA/RData/HARV_SS.RData") # Spectral species
## Spectra data range from column 4 to 430
plotNames <- unique(ss_threshold$plotID)
nPlots <- length(plotNames)
Range <- 5:430 ## Spectra data range from column 5 to 430

x <- demon_momentsSPEC(spectra = ss_threshold, nPlots = nPlots, 
                       plotNames = plotNames, specRange = Range)

load("DATA/Harvard_NEON/Data/harv_matchData.RData") # HARV data
rm(harv_comm, harv_taxonomy)

traits <- matchTrait[, c(1, 10)]
names(traits) <- c("Species", "WPH")

xx <- demon_momentsTraitPhylo(comm = matchPhylo$comm, 
                              trait = traits, 
                              phylo = matchPhylo$phy, 
                              plotNames = plotNames, 
                              nPlots = nPlots)
xx$moments_taxo
xx$moments_traits
xx$moments_phylo
