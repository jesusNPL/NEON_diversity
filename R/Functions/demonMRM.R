
library(tidyverse)


##### Make raster using spectra data #####

# spectra = data.frame with spectra data
# plotNames = plot names
# nPlots = mumber of plots per site

demon_SPEC_to_RASTER <- function(spectra, plotNames, nPlots) { 
  
  if ( ! ("raster" %in% installed.packages())) {install.packages("raster", dependencies = TRUE)} 
  if ( ! ("svMisc" %in% installed.packages())) {install.packages("svMisc", dependencies = TRUE)} 
  
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

multiMRM <- function(spectra, site, distance = "euclidean", iters = 1000) { 
  ## Inits
  plotNames <- na.omit(unique(spectra$plotID))
  
  cbn <- data.frame(apply(combn(plotNames, 2), 2, paste, collapse = '-'))
  names(cbn) <- "comb"
  
  cbn <- cbn %>% 
    mutate(combination = comb) %>% 
    separate(comb, into = c("plot1", "plot2"), sep = "-")
  
  ## Store results
  coefficients <- list()
  estimates <- list()
  F_test <- list()
  
  for(i in 1:nrow(cbn)) { 
    
    print(paste0("Performing pairwise MRM for plots ", cbn[i, 3], " at ", site ," ..."))
    
    p1 <- cbn[i, 1]
    p2 <- cbn[i, 2]
    
    plot1 <- spectra %>% 
      filter(plotID == p1)
    
    plot2 <- spectra %>% 
      filter(plotID == p2)
    
    ## Distances 
    p1_dist <- factoextra::get_dist(plot1[, 5:430], method = distance, stand = TRUE)
    
    p2_dist <- factoextra::get_dist(plot2[, 5:430], method = distance, stand = TRUE)
    
    ### MRM ### 
    fit <- ecodist::MRM(p1_dist ~ p2_dist, nperm = iters)
    
    ## Coefficientes
    coefs <- data.frame(fit$coef) 
    coefs$param <- rownames(coefs)
    coefs$plotID <- cbn[i, 3]
    coefs$Site <- site 
    rownames(coefs) <-  NULL
    coefficients[[i]] <- coefs
    
    ## Estimates R2
    r2 <- data.frame(fit$r.squared) 
    r2$param <- rownames(r2)
    r2$plotID <- cbn[i, 3]
    r2$Site <- site 
    rownames(r2) <- NULL
    estimates[[i]] <- r2
    
    ## F test 
    ftest <- data.frame(fit$F.test) 
    ftest$param <- rownames(ftest)
    ftest$plotID <- cbn[i, 3]
    ftest$Site <- site 
    rownames(ftest) <- NULL 
    F_test[[i]] <- ftest
    
  } 
  
  coefficients <- do.call(rbind, coefficients)
  estimates <- do.call(rbind, estimates)
  F_test <- do.call(rbind,F_test) 
  
  ## save results
  res <- list(coefficients = coefficients, 
              estimates = estimates, 
              F_test = F_test)
  
  return(res)
  
}


multiMRMdist <- function(distances, plotNames, site, iters = 1000) { 
  ## Inits

  cbn <- data.frame(apply(combn(plotNames, 2), 2, paste, collapse = '-'))
  names(cbn) <- "comb"
  
  cbn <- cbn %>% 
    mutate(combination = comb) %>% 
    separate(comb, into = c("plot1", "plot2"), sep = "-")
  
  ## Store results
  coefficients <- list()
  estimates <- list()
  F_test <- list()
  
  for(i in 1:nrow(cbn)) { 
    
    print(paste0("Performing pairwise MRM for plots ", cbn[i, 3], " at ", site ," ..."))
    
    p1 <- cbn[i, 1]
    p2 <- cbn[i, 2]
    
    plot1 <- distances[[p1]] 
    
    plot2 <- distances[[p2]] 
    
    if(nrow(plot2) != nrow(plot1)) { 
      next
      }

    ## Distances 
    p1_dist <- as.dist(plot1)
    
    p2_dist <- as.dist(plot2)
    
    ### MRM ### 
    fit <- ecodist::MRM(p1_dist ~ p2_dist, nperm = iters)
    
    ## Coefficientes
    coefs <- data.frame(fit$coef) 
    coefs$param <- rownames(coefs)
    coefs$plotID <- cbn[i, 3]
    coefs$Site <- site 
    rownames(coefs) <-  NULL
    coefficients[[i]] <- coefs
    
    ## Estimates R2
    r2 <- data.frame(fit$r.squared) 
    r2$param <- rownames(r2)
    r2$plotID <- cbn[i, 3]
    r2$Site <- site 
    rownames(r2) <- NULL
    estimates[[i]] <- r2
    
    ## F test 
    ftest <- data.frame(fit$F.test) 
    ftest$param <- rownames(ftest)
    ftest$plotID <- cbn[i, 3]
    ftest$Site <- site 
    rownames(ftest) <- NULL 
    F_test[[i]] <- ftest
    
  } 
  
  coefficients <- do.call(rbind, coefficients)
  estimates <- do.call(rbind, estimates)
  F_test <- do.call(rbind,F_test) 
  
  ## save results
  res <- list(coefficients = coefficients, 
              estimates = estimates, 
              F_test = F_test)
  
  return(res)
  
}

#x <- multiMRM(spectra = td, site = "TOOL", distance = "euclidean", iters = 100)

#td <- read.csv("Dropbox/Macrosystems_NEON/DATA/Spectra/TOOL_spectra_grid.csv")

#toolPlots <- na.omit(unique(td$plotID))
#nPlots <- length(toolPlots)

#tool_ras <- demon_SPEC_to_RASTER(spectra = td, 
 #                                plotNames = toolPlots, 
  #                               nPlots = nPlots)

#tool_dist <- distSAM(spectraRAS = tool_ras, nPlots = nPlots, plotNames = tooPlots)                               

#tool_mrm <- multiMRMdist(distances = tool_dist, 
 #                        plotNames = toolPlots, 
  #                       site = "TOOL", 
   #                      iters = 1000)

