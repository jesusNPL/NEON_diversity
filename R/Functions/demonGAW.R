
# method = "optimized" or "analytic"
demon_gawdis <- function(spectra, specRange, plotNames, nPlots, 
                         standardize = FALSE, Ncores, method) {
  if ( ! ("gawdis" %in% installed.packages())) {install.packages("gawdis", dependencies = TRUE)} 
  if ( ! ("dplyr" %in% installed.packages())) {install.packages("dplyr", dependencies = TRUE)}
  if ( ! ("tidyr" %in% installed.packages())) {install.packages("tidyr", dependencies = TRUE)} 
  
  require(dplyr)
  require(tidyr)
  # Start Parallel 
  
  spectra <- na.omit(spectra)
  
  # Start parallel computation
  cl <- parallel::makeCluster(Ncores)
  
  gaw <- list()
  
  for(i in 1:nPlots) { 
    
    print(plotNames[i])
    
    spec <- spectra %>% filter(plotID == plotNames[i]) %>% 
      drop_na() # remove NAs 
    
    specSel <- spec[, specRange]
    # Standardizing
    if(standardize == TRUE)
    {
      specSel <- vegan::decostand(specSel, method = "total", MARGIN = 1)
    } else { 
      specSel <- specSel
    }
    
    if(method == "optimized")
    {
      gaw[[i]] <- gawdis::gawdis(specSel, w.type = "optimized", groups.weight = FALSE, silent = TRUE, opti.maxiter = 200) 
    } else { 
      gaw[[i]] <- gawdis::gawdis(specSel, w.type = "analytic", groups.weight = FALSE, silent = TRUE) 
    }
    
    
    print(paste0("GawDist computed for the plot ", plotNames[i], " starting next plot..."))
  } 
  #Stop parallel computation
  parallel::stopCluster(cl = cl)
  
  names(gaw) <- plotNames
  return(gaw)
}
