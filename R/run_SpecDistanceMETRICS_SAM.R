##### Run alpha spectral diversity #####
library(tidyverse)
library(picante)
library(fundiversity)
library(phytools)

source("R/NEON_diversity/R/Functions/AlphaSpecDIV_SAM.R")

### Load data
load("DATA/matchPhyloComm/matchPhyloComm.RData")

### Inits
sites <- names(matched_PhyComm)
nSites <- length(sites)

Range <- 5:430 ## Spectra data range from column 5 to 430
QS <- 0 # different q values

##### Create rasters ######

NEON_ras <- list()

for(j in 1:nSites) {
  print(sites[j])
  
  spec <- read.csv(paste0("DATA/Spectra/", sites[j], "_spectra_grid.csv"))
  
  spec <- spec %>% 
    drop_na() 
  
  plotNames <- unique(spec$plotID)
  nPlots <- length(plotNames)
  
  NEON_ras[[j]] <- demon_SPEC_to_RASTER(spectra = spec, 
                                   plotNames = plotNames, 
                                   nPlots = nPlots)
  
}

names(NEON_ras) <- sites
NEON_ras$ABBY$P_1

##### Make SAM #####

NEON_sam <- list()

for(j in 1:nSites) {
  print(sites[j])
  
  spec <- NEON_ras[[j]]
 
  
  plotNames <- names(spec)
  nPlots <- length(plotNames)
  
  NEON_sam[[j]] <- distSAM(spectraRAS = spec, 
                      plotNames = plotNames, 
                      nPlots = nPlots)
  
}

names(NEON_sam) <- sites
NEON_sam$ABBY$P_1

##### Run spectral diversity based distances ######

NEON_specDiv_SAM <- list()

for(j in 1:nSites) {
  
  specSAM <- NEON_sam[[j]]
  
  plotNames <- names(specSAM)
  nPlots <- length(plotNames)
  
  spec_div <- demon_DivDistance_SAM(samList = specSAM, 
                                    Q = 0, nPlots = nPlots, 
                                    plotNames = plotNames)
  spec_div$Site <- sites[j] 
  
  NEON_specDiv_SAM[[j]] <- spec_div
}

NEON_specDiv_SAM <- do.call(rbind, NEON_specDiv_SAM)

saveRDS(NEON_specDiv_SAM, file = "Results/spectralDiversity/Reanalyses/spectral_NEON_SAM_q0.rds")
