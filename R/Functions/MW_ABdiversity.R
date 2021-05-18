
##### Function to estimate spectral beta diversity and its components #####
# spectra = spectra data
# coords = coordinates for each cell/pixel
# window = window size, e.g., 1.5
# metric = dissimilarity index, it can be "bray" or "ruzicka"

demon_specBeta <- function(spectra, coords, specRange, window, metric, standardize = FALSE, plotNames) { 
  
  spectra <- na.omit(spectra)
  
  spec <- spectra[, specRange]
  
  if(standardize == TRUE)
  {
    spec <- vegan::decostand(spec, method = "total", MARGIN = 1)
  } else { 
    spec <- spec 
  }
  
  coords <- spectra[, coords]
  
  spec <- cbind(coords, spec)
  
  nCells <- nrow(spec)
  
  
  mean_turnover <- numeric(length(spectra[, 1]))
  mean_nestedness <- numeric(length(spectra[, 1]))
  mean_beta <- numeric(length(spectra[, 1]))
  
  for(i in 1:nCells){
    svMisc::progress(i, max.value = nCells)
    
    specSel <- CommEcol::select.window(xf = spec[i, 1], 
                          yf = spec[i, 2], radius = window, 
                          xydata = spec)[, -c(1, 2)]
    
    res <- betapart::beta.pair.abund(specSel, index.family = metric)
    
    mean_turnover[i] <- mean(as.matrix(res[[1]])[2:length(as.matrix(res[[1]])[, 1]), 1], na.rm = TRUE)
    mean_nestedness[i] <- mean(as.matrix(res[[2]])[2:length(as.matrix(res[[2]])[, 1]), 1], na.rm = TRUE)
    mean_beta[i] <- mean(as.matrix(res[[3]])[2:length(as.matrix(res[[3]])[, 1]), 1], na.rm = TRUE)
  }
  # beta data.frame
  specBeta_df <- data.frame(coords, mean_turnover, mean_nestedness, mean_beta)
  specBeta_df$betaRatio <- (mean_nestedness/mean_turnover)
  specBeta_df$plotID <- plotNames
  names(specBeta_df) <- c("X", "Y", "specTurnover", "specNestedness", "specBeta", "specBetaRatio", "plotID")
  # beta spatial
  specBeta_spt <- sp::SpatialPixelsDataFrame(points = specBeta_df[, 1:2], specBeta_df[, 3:6], 
                                             proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  # beta raster
  specBeta_ras <- raster::brick(specBeta_spt) 
  
  results <- list("specBeta" = specBeta_df, 
                  "specBetaPoints" = specBeta_spt, 
                  "specBetaRaster" = specBeta_ras)
  
  return(results) 
}


##### Function to estimate spectral RAO, RICHNESS, EVENNESS and DIVERGENCE ######
# spectra = spectra data
# coords = coordinates for each cell/pixel
# window = window size, e.g., 1.5
# nPCA = number of PCA axes

demon_specFour <- function(spectra, coords, specRange, window, nPCA, standardize = FALSE, plotNames) {
  
  spectra <- na.omit(spectra)
  
  spec <- spectra[, specRange]
  
  if(standardize == TRUE)
  {
    spec <- vegan::decostand(spec, method = "total", MARGIN = 1)
  } else { 
    spec <- spec 
  }
  
  coords <- spectra[, coords]
  
  nCells <- nrow(spec)
  
  spec.pca <- prcomp(spec, center = TRUE, scale. = TRUE)
  
  specPCA_sel <- cbind(coords, data.frame(spec.pca$x[, 1:nPCA]))
  
  SSQ <- numeric(length = nCells)
  SSRic <- numeric(length = nCells)
  SSEve <- numeric(length = nCells)
  SSDiv <- numeric(length = nCells)
  
  for(i in 1:nCells){
    svMisc::progress(i, max.value = nCells)
    
    spec <- CommEcol::select.window(xf = specPCA_sel[i, 1], 
                                    yf = specPCA_sel[i, 2], radius = window, 
                                    xydata = specPCA_sel)[, -c(1, 2)]
    
    SSQ[i] <- fundiversity::fd_raoq(spec)[2]
    SSRic[i] <- fundiversity::fd_fric(spec)[2]
    SSEve[i] <- fundiversity::fd_feve(spec)[2]
    SSDiv[i] <- fundiversity::fd_fdiv(spec)[2]
    
  }
  
  # Diversity measurements
  specDIV_df <- data.frame(coords, unlist(SSQ), unlist(SSRic), unlist(SSDiv), unlist(SSEve))
  specDIV_df$plotID <- plotNames
  names(specDIV_df) <- c("X", "Y", "specRAOq", "specRichness", "specDivergence", "specEvenness", "plotID")
  # Diversity measurements spatial
  specDIV_spt <- sp::SpatialPixelsDataFrame(points = specDIV_df[, 1:2], specDIV_df[, 3:6], 
                                             proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  # Diversity measurements raster
  specDIV_ras <- raster::brick(specDIV_spt) 
  
  # PCA data
  specPCA <- cbind(coords, data.frame(spec.pca$x))
  
  specPCA_spt <- sp::SpatialPixelsDataFrame(points = specPCA[, 1:2], data.frame(spec.pca$x), 
                                            proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  specPCA_ras <- raster::brick(specPCA_spt) 
  
  results <- list("specDiversity" = specDIV_df, 
                  "specDiversityPoints" = specDIV_spt, 
                  "specDiversityRaster" = specDIV_ras, 
                  "specPCA" = specPCA, 
                  "specPCARaster" = specPCA_ras)
  
  return(results) 

}


##### Function to estimate spectral SCHEINER metrics and mapping ######
# spectra = spectra data
# coords = coordinates for each cell/pixel
## specRange = position of the bands in the data frame.
# window = window size, e.g., 1.5
## standardize = logical (TRUE/FALSE). If TRUE pixel values are standardized using the function decostand using the {vegan} package
## gowDist = logical (TRUE/FALSE). If TRUE the distance using GOWER distance is applied.
## distance = the distance measure to be used. IF gowDist = TRUE, distance = NULL
# This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", 
# "minkowski", "pearson", "spearman" or "kendall".
## Q = numeric value that represents the Hill coefficient q
## plotNames = names of the plots


demon_Scheiner_MW <- function(spectra, coords, specRange, window, standardize = FALSE, 
                             distance = "euclidean", gowDist = FALSE, Q, plotNames) {
  library(dplyr)
  library(tidyr)
  library(raster)
  # Scheiner metrics are calculated using functions written by Shan Kothari.
  source("https://raw.githubusercontent.com/ShanKothari/DecomposingFD/master/R/AlphaFD.R")
  
  spectra <- na.omit(spectra)
  
  spec <- spectra[, specRange]
  
  if(standardize == TRUE)
  {
    spec <- vegan::decostand(spec, method = "total", MARGIN = 1)
  } else { 
    spec <- spec 
  }
  
  coords <- spectra[, coords]
  
  spec <- cbind(coords, spec)
  
  nCells <- nrow(spec)
  
  ## getting qDT, qDTM, and qEt from qHt
  nSP <- numeric(length = nCells)
  q <- numeric(length = nCells)
  M <- numeric(length = nCells)
  M_prime <- numeric(length = nCells)
  qHt <- numeric(length = nCells)
  qEt <- numeric(length = nCells)
  qDT <- numeric(length = nCells)
  qDTM <- numeric(length = nCells)
  
  for(i in 1:nCells){
    svMisc::progress(i, max.value = nCells)
    
    specSel <- CommEcol::select.window(xf = spec[i, 1], 
                                    yf = spec[i, 2], radius = window, 
                                    xydata = spec)[, -c(1, 2)]
    # if we wanted to transform our original matrix into relative abundance we
    # can use the same function and the "total" method, which divides each value in a row
    # by the total of all values in a row if the MARGIN is set to ONE.

    if (gowDist == TRUE) {
      specDist <- FD::gowdis(specSel)
    } else {
      specDist <- factoextra::get_dist(specSel, method = distance, stand = TRUE)
    }
    
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
    
  }
  
  # Diversity measurements dataframe
  Scheiner_df <- data.frame(coords, nSP, q, M, M_prime, qHt, qEt, qDT, qDTM)
  Scheiner_df$PlotID <- plotNames
  names(Scheiner_df) <- c("X", "Y", "nPixels", "Q", "M", "M_prime", "qHt", "qEt", "qDT", "qDTM", "plotID")
  # Diversity measurements spatialpoints
  Scheiner_spt <- sp::SpatialPixelsDataFrame(points = Scheiner_df[, 1:2], Scheiner_df[, 5:10], 
                                            proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  # Diversity measurements rasters
  Scheiner_ras <- raster::brick(Scheiner_spt) 
  
  results <- list("specScheiner" = Scheiner_df, 
                  "specScheinerPoints" = Scheiner_spt, 
                  "specScheinerRaster" = Scheiner_ras)
  
  return(results) 

}




