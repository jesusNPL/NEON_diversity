##### Run alpha spectral diversity #####
library(tidyverse)
library(picante)
library(fundiversity)
library(phytools)

source("R/NEON_diversity/R/Functions/AlphaSpecDIV_SAM.R")

### Load data
load("DATA/matchPhyloComm/matchPhyloComm.RData")

### Inits
# Site names
sites <- names(matched_PhyComm)
# Number of plots
nSites <- length(sites)

##### Create rasters ######

NEON_ras <- list()

for(j in 1:nSites) { 
  
  print(sites[j])
  
  # spectra data normalized
  spec <- read.csv(paste0("DATA/Spectra/normalized/normalized_mask_", sites[j], "_spectra_grid.csv"))
  
 # spec <- spec %>% # remove NA pixel values - those removed by the NDVI and NIR masking process
  #  drop_na() 

  # plot names
  plotNames <- unique(spec$plotID)
  
  # number of plots within the NEON site
  nPlots <- length(plotNames)
  
  NEON_ras[[j]] <- demon_SPEC_to_RASTER(spectra = spec, 
                                   plotNames = plotNames, 
                                   nPlots = nPlots)
  
}

names(NEON_ras) <- sites

plot(NEON_ras$OAES$P_1) # plot several bands

plot(NEON_ras$OAES$P_1$X737) # plot specific bands

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

### Save SAM distances 
save(NEON_sam, 
     file = "DATA/Spectral_distance_SAM/NEON_spec_SAM.RData")

##### Run spectral diversity based distances ######

NEON_specDiv_SAM <- list()

for(j in 1:nSites) {
  
  print(sites[j])
  
  specSAM <- NEON_sam[[j]]
  
  plotNames <- names(specSAM) 
  
  nPlots <- length(plotNames)
  
  spec_div <- demon_DivDistance_SAM(samList = specSAM, 
                                    Q = 0, 
                                    nPlots = nPlots, 
                                    plotNames = plotNames)
  spec_div$Site <- sites[j] 
  
  NEON_specDiv_SAM[[j]] <- spec_div 
  
}

### Combine results in a single data.frame
NEON_specDiv_SAM <- do.call(rbind, NEON_specDiv_SAM) 

NEON_specDiv_SAM <- NEON_specDiv_SAM %>% 
  dplyr::select(Site, everything())

### save results
saveRDS(NEON_specDiv_SAM, 
        file = "Results/spectralDiversity/Reanalyses/spectral_normalized_NEON_SAM_q0.rds")

write_csv(NEON_specDiv_SAM, 
          file = "Results/spectralDiversity/Reanalyses/spectral_normalized_NEON_SAM_q0.csv")

