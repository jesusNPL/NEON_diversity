
##### Function to set the number of spectral species #####

# spectra = spectra data
# distance = distance metric, e.g., euclidean
# Nspp = observed number of species
# metrics = vector of metrics to be evaluated
# Ncores = number of cores to be used

demon_ThreshClusters <- function(spectra, distance = "euclidean", Nspp, metrics, Ncores, progress = FALSE) { 
  
  if ( ! ("NbClust" %in% installed.packages())) {install.packages("NbClust", dependencies = TRUE)} 
  if ( ! ("dplyr" %in% installed.packages())) {install.packages("dplyr", dependencies = TRUE)}
  if ( ! ("tidyr" %in% installed.packages())) {install.packages("tidyr", dependencies = TRUE)}
  if ( ! ("svMisc" %in% installed.packages())) {install.packages("svMisc", dependencies = TRUE)} 
  if ( ! ("parallel" %in% installed.packages())) {install.packages("parallel", dependencies = TRUE)} 
  
  require(dplyr)
  require(tidyr)
  # remove NAs
  spec <- na.omit(spectra)
  
  #  Distance
  Dist_mat <- dist(spec, method = distance, diag = FALSE) 
  
  resu <- data.frame(matrix(ncol = 2, nrow = length(metrics)))
  
  # Start Parallel 
  cl <- parallel::makeCluster(Ncores)
  #register it to be used by %dopar%
  #doParallel::registerDoParallel(cl = cl)
  
  for(j in 1:length(metrics)) {
    
    if (progress == TRUE) { 
      
      svMisc::progress(j, max.value = length(metrics))

    }
    
    res <- NbClust::NbClust(data = spec, 
                            diss = Dist_mat, 
                            distance = NULL, 
                            min.nc = 3, max.nc = Nspp, 
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
  
  resu$metric <- metrics
  names(resu) <- c("Number_clusters", "Value_Index", "Metric") 

  # Majority rule
  mrClust <- resu %>% 
    count(Number_clusters) #%>% 
  #filter(n > 1)
  meanClust <- round(mean(resu$Number_clusters)) 
  
  results <- list("Metrics" = resu, "Majority_rule" = mrClust, "Mean_rule" = meanClust)
  
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
    
    if (threshold[i] <= 3) {
      threshold[i] <- 4
    }
    
    th <- threshold[i]
    
    specPlot <- spectra %>% filter(plotID == plotNames[i]) 
    # remove NAs
    spec <- na.omit(specPlot[, 5:430])
    # Scale spectra
    spec <- scale(spec) 
    # Perform  robust version of K-means
    robustK <- pam(specPlot, th) 
    
    specPlot$SS <- as.factor(robustK$clustering)
    
    SS[[i]] <- specPlot
    
  }
  
  SS <- data.frame(do.call(rbind, SS))
  
  return(SS)
  
}
