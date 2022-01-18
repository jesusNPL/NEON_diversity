##### Make raster using spectra data #####

# spectra = data.frame with spectra data
# plotNames = plot names
# nPlots = mumber of plots per site

demon_SPEC_to_RASTER <- function(spectra, plotNames, nPlots) { 
  
  if ( ! ("raster" %in% installed.packages())) {install.packages("raster", dependencies = TRUE)} 
  if ( ! ("dplyr" %in% installed.packages())) {install.packages("dplyr", dependencies = TRUE)}
  if ( ! ("tidyr" %in% installed.packages())) {install.packages("tidyr", dependencies = TRUE)}
  if ( ! ("svMisc" %in% installed.packages())) {install.packages("svMisc", dependencies = TRUE)} 
  
  require(dplyr)
  require(tidyr)
  require(raster)
  
  plots <- list()
  
  for(i in 1:nPlots) { 
    svMisc::progress(i, max.value = nPlots)
    
    pt <- spectra %>% filter(plotID == plotNames[i])
    pt_spt <- sp::SpatialPixelsDataFrame(points = pt[, 1:2], pt[, 5:430])
    
    tmp_ras <- raster::brick(pt_spt)
    plots[[i]] <- tmp_ras
    
  }
  names(plots) <- plotNames
  return(plots)
}

#harv_ras <- demon_SPEC_to_RASTER(spectra = harv, 
 #                                plotNames = harvPlots, 
  #                               nPlots = nPlots)


##### calculate SAM ######
distSAM <- function(spectraRAS, nPlots, plotNames) {
  
  library(RStoolbox)

  sam_lst <- list()
  
  for(i in 1:nPlots) { 
    print(plotNames[i])
    
    spec <- spectraRAS[[i]]
    spec_coords <- sp::coordinates(spec) 
    endmembers <- raster::extract(spec, spec_coords)
    sam_Ras <- RStoolbox::sam(spec, endmembers, angles = TRUE)
    sam_DT <- as.data.frame(sam_Ras)
    sam_lst[[i]] <- sam_DT
  }
  names(sam_lst) <- plotNames
  return(sam_lst)
}


##### Function to estimate spectral diversity using distance metrics based SAM #####
## spectra = spectra data
## standardize = logical (TRUE/FALSE). If TRUE pixel values are standardized using the function decostand using the {vegan} package
## gowDist = logical (TRUE/FALSE). If TRUE the distance using GOWER distance is applied.
## distance = the distance measure to be used. IF gowDist = TRUE, distance = NULL
# This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", 
# "minkowski", "pearson", "spearman" or "kendall".
## Q = numeric value that represents the Hill coefficient q
## nPlots = number of plots
## plotNames = names of the plots
## specRange = position of the bands in the data frame.

demon_DivDistance_SAM <- function(samList, Q, nPlots, plotNames) {
  library(dplyr)
  library(tidyr)
  # Scheiner metrics are calculated using functions written by Shan Kothari.
  source("https://raw.githubusercontent.com/ShanKothari/DecomposingFD/master/R/AlphaFD.R")
  
  ### Store calculations ###
  SD <- numeric(length = nPlots) # analogous to PD - phylogenetic diversity
  MSD <- numeric(length = nPlots) # analogous to MPD - mean pairwise distance
  MS <- numeric(length = nPlots) # simple mean distance
  MNSD <- numeric(length = nPlots) # analogous to MNTD - mean nearest taxon distance
  MxNSD <- numeric(length = nPlots) # analogous to MNTD - mean nearest taxon distance
  SNTD <- numeric(length = nPlots) # analogous to SNTD - standard deviation Of MNTD
  aveRange <- numeric(length = nPlots) # Average range
  
  ## getting qDT, qDTM, and qEt from qHt
  nSP <- numeric(length = nPlots)
  q <- numeric(length = nPlots)
  M <- numeric(length = nPlots)
  M_prime <- numeric(length = nPlots)
  qHt <- numeric(length = nPlots)
  qEt <- numeric(length = nPlots)
  qDT <- numeric(length = nPlots)
  qDTM <- numeric(length = nPlots)
  
  # Functional Divergence (FDiv)
  SDiv <- numeric(length = nPlots)
  # Functional Evenness (FEve)
  SEve <- numeric(length = nPlots)
  # Rao's entropy index (Rao's Q)
  SRaoQ <- numeric(length = nPlots)
  # Functional Dispersion (FDis)
  SDis <- numeric(length = nPlots)
  
  for(i in 1:nPlots) { 
    
    svMisc::progress(i, max.value = nPlots)
    
    spec <- samList[[i]]
    
    spec <- na.omit(spec) 
    
    if(nrow(spec) == 0) {
      next
    }
    
    specDist <- as.dist(spec)

    # if we wanted to transform our original matrix into relative abundance we
    # can use the same function and the "total" method, which divides each value in a row
    # by the total of all values in a row if the MARGIN is set to ONE.
    
    ### Start calculations ###
    SD[i] <- sum(specDist) # analogous to PD - phylogenetic diversity
    msd <- specDist[lower.tri(specDist, diag = FALSE)]
    MSD[i] <- mean(na.omit(msd)) # analogous to MPD - mean pairwise distance
    MS[i] <- mean(na.omit(specDist)) # simple mean distance
    #mnsd <- specDist
    #diag(mnsd) <- NA
    #MNSD[i] <- mean(apply(mnsd, MARGIN = 2, min), na.rm = TRUE) # analogous to MNTD - mean nearest taxon distance
    #MxNSD[i] <- max(apply(mnsd, MARGIN = 2, min), na.rm = TRUE) # maximum MNTD
    #SNTD[i] <- sd(apply(mnsd, MARGIN = 2, min), na.rm = TRUE) # standard deviation Of MNTD
    aveRange[i] <- mean(apply(spec, MARGIN = 2, max) - apply(spec, MARGIN = 2, min))
    
    # Scheiner metrics
    Scheiner <- FTD(specDist, q = Q)
    nSP[i] <- Scheiner[[1]]
    q[i] <- Scheiner[[2]]
    M[i] <- Scheiner[[3]]
    M_prime[i] <- Scheiner[[4]]
    qHt[i] <- Scheiner[[5]]
    qEt[i] <- Scheiner[[6]]
    qDT[i] <- Scheiner[[7]]
    qDTM[i] <- Scheiner[[8]]
    
    # Compute Functional Divergence (FDiv)
    SDiv[i] <- fundiversity::fd_fdiv(spec)[2]
    # Compute Functional Evenness (FEve)
    SEve[i] <- fundiversity::fd_feve(spec)[2]
    # Compute Rao's entropy index (Rao's Q)
    SRaoQ[i] <- fundiversity::fd_raoq(spec)[2]
    # Compute Functional Dispersion (FDis)
    SDis[i] <- fundiversity::fd_fdis(spec)[2]
    
  }
  
  alphaSDiv <- data.frame(plotNames, SD, MSD, MS, MNSD, MxNSD, SNTD, aveRange,  
                          nSP, q, M, M_prime, qHt, qEt, qDT, qDTM, 
                          unlist(SDiv), unlist(SEve), unlist(SRaoQ), unlist(SDis)) 
  
  names(alphaSDiv) <- c("plotID", "SD", "MSD", "MSm", "MNSD", "MxNSD", "SNTD", "Range", 
                        "nPixels", "qHill", "M", "mPrime", "qHt", "qEt", "qDT", "qDTM", 
                        "SDivergence", "SEvenness", "SRao", "SDispersion")
  
  return(alphaSDiv)
  
}

