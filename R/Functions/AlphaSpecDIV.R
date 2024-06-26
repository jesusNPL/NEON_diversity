
#load("DATA/RData/HARV_SS.RData")

#harv_p1 <- ss_NO_threshold %>% 
 # filter(plotID == "P_1")

##### Function to estimate spectral diversity using distance metrics #####
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

demon_DivDistance <- function(spectra, distance, standardize = FALSE, gowDist = FALSE, 
                              Q, nPlots, plotNames, specRange) {
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
    
    spec <- spectra %>% filter(plotID == plotNames[i]) %>% 
      na.omit()
    
    spec <- spec[, specRange]
    
    # if we wanted to transform our original matrix into relative abundance we
    # can use the same function and the "total" method, which divides each value in a row
    # by the total of all values in a row if the MARGIN is set to ONE.
    if(standardize == TRUE)
    {
      spec <- vegan::decostand(spec, method = "total", MARGIN = 1)
    } 
    
    if (gowDist == TRUE) {
      specDist <- FD::gowdis(spec)
    } else {
      specDist <- factoextra::get_dist(spec, method = distance, stand = TRUE)
    }
    
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

##### Function to calculate phylogenetic diversity ##### 

## comm: picante community data matrix (CDM)
## phylo: phylogenetic tree that match the names of the species in the CDM
## Q = numeric value that represents the Hill coefficient q
## plotNames = names of the plots

demon_PhyloDistance <- function(comm, phylo, abundance = TRUE, 
                                Q, siteName, plotNames) {
  # Scheiner metrics are calculated using functions written by Shan Kothari.
  source("https://raw.githubusercontent.com/ShanKothari/DecomposingFD/master/R/AlphaFD.R")
  library(picante)
  
  ### Scheiner metrics 
  Scheiner <- FTD.comm(tdmat = cophenetic(phylo), 
                       spmat = comm, 
                       q = Q, 
                       abund = TRUE)$com.FTD 
  
  print("Scheiner metrics, done!")

  ### Picante metrics
  PD <- picante::ses.pd(samp = comm, 
                        tree = phylo, 
                        null.model = "taxa.labels")  
  
  MPD <- picante::ses.mpd(samp = comm, cophenetic(phylo), 
                          null.model = "taxa.labels", 
                          abundance.weighted = TRUE) 
  
  PSDs <- psd(samp = comm, 
              tree = phylo)
    
  print("Picante metrics, done!")
  
  # PD and MPD
  PICANTE <- data.frame(PD[, c("pd.obs", "pd.obs.z")], 
                        MPD[, c("mpd.obs", "mpd.obs.z")], 
                        PSDs[, c("PSV", "PSC", "PSR", "PSE")]) 
  
  ### Rao Q
  PICANTE$RaoD <- raoD(comm, force.ultrametric(phylo))$Dkk
  
  print("Rao Q diversity metric, done! ")

  ### Store results
  res <- data.frame(siteID = siteName, 
                    plotID = plotNames, 
                    PICANTE, 
                    Scheiner)
  
  names(res) <- c("siteID", "plotID", "PD", "PDz", "MPD", "MPDz", 
                  "PSV", "PSC", "PSR", "PSE", "RaoD", 
                  "SR", "qHill", "M", "mPrime", "qHt", "qEt", "qDT", "qDTM")
  
  return(res)

}

##### Function to calculate trait diversity ##### 

## comm: picante community data matrix (CDM)
## trait: traits that match the names of the species in the CDM
## Q = numeric value that represents the Hill coefficient q
## plotNames = names of the plots
## gowDist = logical (TRUE/FALSE). If TRUE the distance using GOWER distance is applied.

demon_TraitDistance <- function(comm, trait, abundance = TRUE, Q, plotNames, nPlots, 
                                gowDist = TRUE, gawDist = FALSE) {
  # Scheiner metrics are calculated using functions written by Shan Kothari.
  source("https://raw.githubusercontent.com/ShanKothari/DecomposingFD/master/R/AlphaFD.R")
  library(picante)
  library(gawdis)
  
  if (gowDist == TRUE) {
    traitDist <- FD::gowdis(trait)
  } else if (gawDist == TRUE) {
    traitDist <- gawdis(trait, w.type = "optimized") 
  } else {
    traitDist <- factoextra::get_dist(trait, method = "euclidean", stand = TRUE)
  }
  
  #matched <- match.comm.dist(comm = matchPhylo$comm, dis = traitDist)
  
  # Scheiner metrics 
  Scheiner <- FTD.comm(tdmat = traitDist, 
                       spmat = comm, 
                       q = Q, 
                       abund = TRUE, 
                       match.names = TRUE)$com.FTD
  
  print("Scheiner metrics, done!")
  
  ### Picante metrics
  MPD <- picante::ses.mpd(samp = comm, traitDist, 
                          null.model = "taxa.labels", abundance.weighted = TRUE) 
  
  PICANTE <- data.frame(MPD[, c("mpd.obs", "mpd.obs.z")]) 
  
  print("Picante metrics, done!")
  
  ## getting frao, fric, feve, fdiv, fdis
  nSP <- numeric(length = nPlots) # Number of species
  trao <- numeric(length = nPlots) # Compute Rao's entropy index (Rao's Q)
  tric <- numeric(length = nPlots) # Compute Functional Richness (FRic)
  teve <- numeric(length = nPlots) # Compute Functional Evenness (FEve)
  tdiv <- numeric(length = nPlots) # Compute Functional Divergence (FDiv)
  tdis <- numeric(length = nPlots) # Compute Functional Dispersion (FDis) 
  
  # start calculations
  comm2 <- as.matrix(comm)
  
  for(i in 1:nPlots) {
    spNames <- comm2[i, ] > 0 
    selection <- names(spNames[spNames == TRUE])
    
    nSP[i] <- length(selection)
    trait_comm <- trait[rownames(trait) %in% selection, ]
    trao[i] <- as.numeric(fd_raoq(trait_comm)[2])
    tric[i] <- as.numeric(fd_fric(trait_comm)[2])
    teve[i] <- as.numeric(fd_feve(trait_comm)[2])
    tdiv[i] <- as.numeric(fd_fdiv(trait_comm)[2])
    tdis[i] <- as.numeric(fd_fdis(trait_comm)[2])

  }
  
  traitDiver <- data.frame(tdiv, teve, trao, tdis, tric)
  
  print("Classic trait diversity metrics, done! ")
  
  ### Store results
  results <- data.frame(plotID = plotNames, PICANTE, Scheiner, traitDiver)
  
  names(results) <- c("plotID", "MFD", "MFDz", # Picante 
                  "SR", "qHill", "M", "mPrime", "qHt", "qEt", "qDT", "qDTM", # Scheiner
                  "TDivergence", "TEveness", "TRao", "TDispersion", "TRichness") # classic
  
  return(results)
  
}



##### Function to calculate RAO and EVE #####

demon_SpecRaoEve <- function(spectra, standardize = FALSE, 
                             nPlots, plotNames, specRange) {
  
  library(dplyr)
  library(tidyr)
  
  ### Store calculations ##### 
  #Dkk <- NULL # within community diversity
  #Dkl <- NULL # Among-community diversity
  #H <- NULL # Among-community diversities excluding within-community diversity
  
  # From RAO function in PICANTE
  totalSpec <- numeric(length = nPlots) # total diversity
  alphaSpec <- numeric(length = nPlots) # alpha diversity - average within pixels (average Dkk)
  betaSpec <- numeric(length = nPlots) # beta diversity - average within pixels (average Dkl)
  FstSpec <- numeric(length = nPlots) # ratio between Beta diversity/total diversity
  
  # From EVE function in FUNDIV
  Renji <- numeric(length = nPlots) 
  invSimpson <- numeric(length = nPlots)
  Pielou <- numeric(length = nPlots)
  Evar <- numeric(length = nPlots)
  
  for(i in 1:nPlots) { 
    
    svMisc::progress(i, max.value = nPlots)
    
    spec <- spectra %>% filter(plotID == plotNames[i]) %>% 
      na.omit()
  
    spec <- spec[, specRange]
    
    # if we wanted to transform our original matrix into relative abundance we
    # can use the same function and the "total" method, which divides each value in a row
    # by the total of all values in a row if the MARGIN is set to ONE.
    if(standardize == TRUE)
    {
      spec <- vegan::decostand(spec, method = "total", MARGIN = 1)
    } else { 
      spec <- spec
      }
    
    ### Start calculations ###
    
    # RAO {picante}
    picanteRAO <- picante::raoD(spec)
    
    #Dkk <- PicanteRAO$Dkk 
    #Dkl <- apply(PicanteRAO$Dkl, MARGIN = 2, mean) 
    #H <- apply(PicanteRAO$H, MARGIN = 2, mean) 
    
    totalSpec[i] <- picanteRAO$total
    alphaSpec[i] <- picanteRAO$alpha
    betaSpec[i] <- picanteRAO$beta 
    FstSpec[i] <- picanteRAO$Fst 
    
    # EVE {fundiv}
    fundivEVE <- fundiv::Eve(spec)
    
    Renji[i] <- mean(fundivEVE$Ea1, na.rm = TRUE)
    invSimpson[i] <- mean(fundivEVE$EinvD, na.rm = TRUE)
    Pielou[i] <- mean(fundivEVE$EJ, na.rm = TRUE)
    Evar[i] <- mean(fundivEVE$Evar, na.rm = TRUE)

  } 
  
  metrics <- data.frame(plotNames, totalSpec, alphaSpec, betaSpec, FstSpec, 
                        Renji, invSimpson, Pielou, Evar) 
  names(metrics) <- c("plotID", "raoTotal", "raoAlpha", "raoBeta", "raoRatio", 
                      "Renji", "SimpsonD", "EJ", "Evar")
  
  return(metrics)
}


##### Function to estimate spectral RAO and EVE diversity and its components #####
# spectra = spectra data
# coords = coordinates for each cell/pixel
# specRange = position of the bands on the data.frame
# standardize = logical (TRUE/FALSE). If TRUE pixels are standardized using the function decostand from {vegan}

demon_SpecRaoEve_RASTER <- function(spectra, coords, specRange, 
                                    standardize = FALSE, plotName) {
  
  library(vegan)
  library(raster)
  
  spectra <- na.omit(spectra)
  
  spec <- spectra[, specRange]
  coords <- spectra[, coords]
  
  # if we wanted to transform our original matrix into relative abundance we
  # can use the same function and the "total" method, which divides each value in a row
  # by the total of all values in a row if the MARGIN is set to ONE.
  if(standardize == TRUE)
  {
    spec <- vegan::decostand(spec, method = "total", MARGIN = 1)
  } else { 
    spec <- spec
  }
  
  ### Start calculations ###
  
  # RAO {picante}
  picanteRAO <- picante::raoD(spec)
  
  Dkk <- picanteRAO$Dkk 
  Dkl <- picanteRAO$Dkl
  H <- picanteRAO$H
  rao <- data.frame(Dkk, Dkl, H)
  
  # EVE {fundiv}
  fundivEVE <- fundiv::Eve(spec)
  
  print(paste0("RAO and EVE computed for ", plotName, "... creating rasters!"))
  
  # Combine results
  metrics <- cbind(coords, rao, fundivEVE)
  metrics$plotID <- plotName
  names(metrics) <- c("X", "y", "Dkk", "Dkl", "H", "EinvD", "EJ", "Evar", "Ea1", "plotID")
  
  # Transform to spatial points
  metrics_spt <- sp::SpatialPixelsDataFrame(points = metrics[, 1:2], metrics[, 3:9], 
                                            proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  #create raster brick
  metrics_ras <- raster::brick(metrics_spt) 
  
  results <- list("Alpha_df" = metrics, 
                  "Alpha_spt" = metrics_spt, 
                  "Alpha_ras" = metrics_ras)
  return(results)
  
} 
