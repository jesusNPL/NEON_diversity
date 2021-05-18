##### Transform all plots to raster ##### 

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
    
    pt <- harv_spec %>% filter(plotID == plotNames[i])
    pt_spt <- sp::SpatialPixelsDataFrame(points = pt[, 1:2], pt[, 5:430])
    
    tmp_ras <- raster::brick(pt_spt)
    plots[[i]] <- tmp_ras
    
  }
  names(plots) <- plotNames
  return(plots)
}

##### Calculate different metrics of spectral diversity ##### 

# raster_lst = list of rasters
# nPlots = mumber of plots per site
# dists = distance to to used, e.g., "euclidean
# nBands = number bands per plot
# metric = metric of diversity
# NC = number of cores
# window = window size

demon_RAOdiv <- function(raster_lst, nPlots, nBands, dists, window, metric, NC) {
  
  if ( ! ("pbapply" %in% installed.packages())) {install.packages("pbapply", dependencies = TRUE)} 
  if ( ! ("rasterdiv" %in% installed.packages())) {install.packages("rasterdiv", dependencies = TRUE)} 
  if ( ! ("doParallel" %in% installed.packages())) {install.packages("doParallel", dependencies = TRUE)} 
  
  require(raster)
  require(rasterdiv)
  require(doParallel)
  
  rastDiv <- list() 
  
  for(i in 1:nPlots) { 
    svMisc::progress(i, max.value = nPlots)
    
    pt_tm <- raster_lst[[i]]
    
    unstackRts <- raster::unstack(pt_tm) 
    unstackRts <- unstackRts[1:nBands]
    # Parametric Rao's quadratic entropy with alpha ranging from 1 to 5
    #if(metric == "Rao") { 
    rao <- pbapply::pblapply(unstackRts, rasterdiv::paRao, dist_m = dists, window = window, np = NC)
    
    rastDiv[[i]] <- rao
  } 
  return(rastDiv)
}

# raster_lst = list of rasters
# nPlots = mumber of plots per site
# nBands = number bands per plot
# metric = metric of diversity
# NC = number of cores
# window = window size

demon_RASdiv <- function(raster_lst, nPlots, nBands,  metric, window, NC) {
  
  if ( ! ("pbapply" %in% installed.packages())) {install.packages("pbapply", dependencies = TRUE)} 
  if ( ! ("rasterdiv" %in% installed.packages())) {install.packages("rasterdiv", dependencies = TRUE)} 
  if ( ! ("doParallel" %in% installed.packages())) {install.packages("doParallel", dependencies = TRUE)} 
  
  require(raster)
  require(rasterdiv)
  require(doParallel)
  
  rastDiv <- list() 
  
  for(i in 1:nPlots) { 
    svMisc::progress(i, max.value = nPlots)
    
    pt_tm <- raster_lst[[i]]
    
    unstackRts <- raster::unstack(pt_tm) 
    unstackRts <- unstackRts[1:nBands]
    # Parametric Rao's quadratic entropy with alpha ranging from 1 to 5
    #if(metric == "Rao") { 
    # Pielou's Evenness
    if(metric == "Pielou") { 
      piel <- pbapply::pblapply(unstackRts, rasterdiv::Pielou, window = window, np = NC)
      
      rastDiv[[i]] <- brick(piel)
    } 
    # Rényi's Index
    if(metric == "Renyi") { 
      ren <- pbapply::pblapply(unstackRts, rasterdiv::Renyi, window = window, np = NC)
      
      rastDiv[[i]] <- ren
    } 
    # Shannon's Diversity
    if(metric == "Shannon") { 
      sha <- pbapply::pblapply(unstackRts, rasterdiv::Shannon, window = window, np = NC) 
      
      rastDiv[[i]] <- brick(sha)
    } 
    # Berger-Parker's Index
    if(metric == "BergerParker") { 
      BP <- pbapply::pblapply(unstackRts, rasterdiv::BergerParker, window = window, np = NC) 
      
      rastDiv[[i]] <- brick(BP)
    } 
    # Hill's numbers 
    if(metric == "Hill") { 
      hill <- pbapply::pblapply(unstackRts, rasterdiv::Hill, window = window, np = NC) 
      
      rastDiv[[i]] <- hill
    } 
  } 
  return(rastDiv)
}

##### Aggregate metrics #####

# metric = list of list with RAO estimations
# nPlots = mumber of plots per site
# nBands = number bands per plot 
# plotNames = plot names

extractRAO <- function(metric, nPlots, nBands, plotNames) {
 
  plots <- list()
  bands <- list()
  
  for(i in 1:nPlots) { 
    tm <- metric[[i]] 
    for(j in 1:nBands) {
      rrr <- tm[[j]]
      bands[[j]] <- rrr$window.3$alpha.1
    }
    plots[[i]] <- brick(bands)
  }
  names(plots) <- plotNames
  return(plots)
}

#rrrr <- extractRAO(metric = rao, nPlots = 3, nBands = 10, plotNames = pnames[1:3])

# metric = list of raster bricks with RAO estimations
# nPlots = mumber of plots per site
# plotNames = plot names

raoDIV_agg <- function(metric, nPlots, plotNames) { 
  
  RAOdiv_mean <- list()
  RAOdiv_SD <- list()
  RAOdiv_max <- list()
  RAOdiv_min <- list()
  RAOdiv_var <- list() 
  
  for(r in 1:nPlots) {
    tmp <- metric[[r]] 
    
    RAOdiv_mean[[r]] <- calc(tmp, fun =  mean)
    RAOdiv_SD[[r]] <- calc(tmp, fun =  sd)
    RAOdiv_max[[r]] <- calc(tmp, fun =  max)
    RAOdiv_min[[r]] <- calc(tmp, fun =  min)
    RAOdiv_var[[r]] <- (calc(tmp, fun =  sd)/calc(tmp, fun =  mean))
  }
  names(RAOdiv_mean) <- plotNames
  names(RAOdiv_SD) <- plotNames
  names(RAOdiv_max) <- plotNames
  names(RAOdiv_min) <- plotNames
  names(RAOdiv_var) <- plotNames
  
  results <- list(RAO_mean = RAOdiv_mean, RAO_SD = RAOdiv_SD, 
                  RAO_max = RAOdiv_max, RAO_min = RAOdiv_min, RAO_var = RAOdiv_var) 
  return(results)
  
}

#ttt <- raoDIV_agg(metric = rrrr, nPlots = 3, plotNames = pnames[1:3])
