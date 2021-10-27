##### Run alpha spectral diversity #####
library(tidyverse)
library(picante)
library(fundiversity)
library(phytools)

source("R/NEON_diversity/R/Functions/AlphaSpecDIV.R")

### Inits
Range <- 5:430 ## Spectra data range from column 5 to 430
QS <- 0 # different q values
sites <- names(matched_PhyComm)
nSites <- length(sites)

##### Run spectral diversity based distances ######

specDiv_NEON <- list()

for(j in 1:nSites) {
  
  spec <- read.csv(paste0("DATA/Spectra/", sites[j], "_spectra_grid.csv"))
  spec <- spec %>% 
    drop_na() 
  
  plotNames <- unique(spec$plotID)
  nPlots <- length(plotNames)
  
  spec_div <- demon_DivDistance(spectra = spec, plotNames = plotNames, 
                                     nPlots = nPlots, specRange = Range, Q = QS, 
                                     standardize = TRUE, gowDist = TRUE, 
                                     distance = "manhattan") 
  spec_div$Site <- sites[j] 
  
  specDiv_NEON[[j]] <- spec_div
}

specDiv_NEON <- do.call(rbind, specDiv_NEON)

saveRDS(specDiv_NEON, file = "Results/spectralDiversity/spectral_NEON_q0.rds")

##### Run spectral moments #####
source("R/NEON_diversity/R/Functions/demon_moments.R")

specMoments_NEON <- list()

for(i in 1:nSites) {
  #print(i)
  
  spec <- read.csv(paste0("DATA/Spectra/", sites[i], "_spectra_grid.csv"))
  spec <- spec %>% 
    drop_na() 
  
  plotNames <- unique(spec$plotID)
  nPlots <- length(plotNames)
  
  spec_moment <- demon_momentsSPEC(spectra = spec, nPlots = nPlots, 
                                   plotNames = plotNames, specRange = Range)
  spec_moment$Site <- sites[i] 
  
  specMoments_NEON[[i]] <- spec_moment
}
  
specMoments_NEON <- do.call(rbind, specMoments_NEON)
  
saveRDS(specMoments_NEON, file = "Results/spectralDiversity/spectral_NEON_moments.rds")
