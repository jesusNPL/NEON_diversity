
##### Function to set the number of spectral species #####

# spectra = spectra data
# distance = distance metric, e.g., euclidean
# Nspp = observed number of species
# metrics = vector of metrics to be evaluated
# Ncores = number of cores to be used

demon_ThreshClusters <- function(spectra, 
                                 distance = "euclidean", 
                                 Nspp, 
                                 metrics, 
                                 nMetrics, 
                                 Ncores, 
                                 progress = FALSE) { 
  
  if ( ! ("NbClust" %in% installed.packages())) {install.packages("NbClust", dependencies = TRUE)} 
  if ( ! ("dplyr" %in% installed.packages())) {install.packages("dplyr", dependencies = TRUE)}
  if ( ! ("tidyr" %in% installed.packages())) {install.packages("tidyr", dependencies = TRUE)}
  if ( ! ("svMisc" %in% installed.packages())) {install.packages("svMisc", dependencies = TRUE)} 
  if ( ! ("parallel" %in% installed.packages())) {install.packages("parallel", dependencies = TRUE)} 
  
  require(dplyr)
  require(tidyr)
  # remove NAs
  #spec <- na.omit(spectra)
  
 
  #  Distance
  spec <- spectra # objects from the SAM
  
  # Matrix to store results
  resu <- data.frame(matrix(ncol = 2, nrow = nMetrics))
  names(resu) <- c("N_clusters", "Index_value")
  
  # Start Parallel 
  cl <- parallel::makeCluster(Ncores)
  #register it to be used by %dopar%
  #doParallel::registerDoParallel(cl = cl)
  
  for(j in 1:nMetrics) {
    
    if (progress == TRUE) { 
      
      svMisc::progress(j, max.value = length(metrics))

    }
    
    res <- NbClust::NbClust(data = spec, 
                            diss = NULL,  
                            distance = "euclidean", 
                            min.nc = 2, # minimum number of species present in the plot
                            max.nc = Nspp, # maximum number of species present in the plot
                            method = "kmeans", 
                            index = metrics[j], # error with ccc, scott, marriot, trcovw, tracew
                            #friedman, rubin, 
                            # gamma, gplus and tau are time consuming 
    )
    
    #fviz_nbclust(res)
    res <- res$Best.nc
    resu[j, ] <- res
  }
  #Stop parallel computation
  parallel::stopCluster(cl = cl)
  
  resu$Metric <- metrics
  
  results <- resu
  
  return(results)
  
}

##### Set spectral species #####

# spectra = spectra data
# threshold = number of clusters per plot
# nPlots = number of plots
# plotNames = name of each plot

demon_SS <- function(spectra, threshold, nPlots, plotNames, progress = FALSE) { 
  
  if ( ! ("cluster" %in% installed.packages())) {install.packages("cluster", dependencies = TRUE)} 
  if ( ! ("dplyr" %in% installed.packages())) {install.packages("dplyr", dependencies = TRUE)}
  if ( ! ("tidyr" %in% installed.packages())) {install.packages("tidyr", dependencies = TRUE)}
  if ( ! ("svMisc" %in% installed.packages())) {install.packages("svMisc", dependencies = TRUE)} 
  
  SS <- list()
  
  for(i in 1:nPlots) {
    
    if (progress == TRUE) { 
      
      svMisc::progress(j, max.value = length(metrics))
      
    }
    
    ## Threshold
    th <- threshold[i]
    
    ## Select plots
    specPlot <- spectra %>% 
      filter(plotID == plotNames[i]) 
    
    ## remove NAs that correspond to masked pixels  
    specPlot <- specPlot %>% 
      drop_na() 

    # Perform  robust version of K-means
    robustK <- cluster::pam(x = specPlot[, 4:ncol(specPlot)], th) 
    
    ## Assign spectral species names 
    specPlot <- specPlot %>% 
      mutate(SS = paste0("SS_", as.factor(robustK$clustering))) %>% 
      select(X, Y, plotID, SS, everything())
    
    SS[[i]] <- specPlot
    
  }
  ## Combine results by site
  SS <- data.frame(do.call(rbind, SS))
  
  return(SS)
  
}
