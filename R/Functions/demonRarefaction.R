library(tidyverse)
library(iNEXT)

## Function to run Rarefaction/Extrapolation for multiple communities using the
## function iNEXT
demoniNEXT <- function(samples, Q, nBoots) {
  
  sites <- names(samples)
  nSites <- length(sites)
  
  estimators <- list() 
  objects <- list() 
  
  for(i in 1:nSites) {
    
    print(sites[i]) 
    
    ## Select sites 
    tmp <- t(samples[[i]])
    tmp <- as.data.frame(tmp)
    
    ## Prepare data
    occ <- data.frame(apply(X = tmp, MARGIN = 2, FUN = as.integer)) 
    rownames(occ) <- rownames(tmp) 
    
    abun <- colSums(occ)
    
    abun[abun == 0] <- NA
    
    abun <- na.omit(abun)
    colSel <- names(abun) 
    
    ## Remove assemblages with no species
    occ2 <- occ %>% 
      select(all_of(colSel)) 
    
    ## Run iNEXT
    out <- iNEXT(occ2, q = Q, 
                 datatype = "abundance", 
                 se = TRUE, conf = 0.95, nboot = nBoots) 
    
    ## Store results
    objects[[i]] <- out
    
    est <- out$AsyEst
    names(est) <- c("plotID", "Diversity", "Observed", "Estimator", "s.e.", "LCL", "UCL")
    est$Site <- sites[i]
    
    estimators[[i]] <- est 
    
    gc()
    
  } 
  
  estimators <- do.call(rbind, estimators) 
  
  res <- list(Estimators = estimators, objects = objects)
  
  return(res)
}

## example 
#sampleCompletenessNEON <- demoniNEXT(samples = matchedTrait$commNEON, 
 #                                    Q = c(0, 1, 2), 
  #                                   nBoots = 1000)

##### Function to get tables #####
getSampleCoverage <- function(iNEXT_object, sites) { 
  
  names(iNEXT_object) <- sites
  
  nSites <- length(sites)
  
  sampleCover <- list() 
  
  for(i in 1:nSites) { 
    
    obj <- iNEXT_object[[i]]$iNextEst 
    
    obj <- do.call(rbind, obj)
    obj$plot <- rownames(obj)
    
    obj <- obj %>% 
      separate(plot, into = c("plotID", "repetition"), sep = "[.Double]") %>% 
      filter(method == "observed")
    
    obj$Site <- sites[i] 
    rownames(obj) <- NULL
    
    sampleCover[[i]] <- obj
    
  }
  
  sampleCover <- do.call(rbind, sampleCover)
  return(sampleCover)
  
}

## Example
#sampCover_NEON <- getSampleCoverage(iNEXT_object = sampleCompletenessNEON_q0$objects, 
 #                                   sites = sites)

##### Rarefaction using vegan #####
demonRarefy <- function(samples) {
  
  sites <- names(samples)
  nSites <- length(sites)
  
  estimators <- list()
  
  for(i in 1:nSites) {
    
    print(sites[i]) 
    
    ## Select sites 
    tmp <- samples[[i]]
    tmp <- as.data.frame(tmp)
    
    ## Prepare data
    occ <- data.frame(apply(X = tmp, MARGIN = 2, FUN = as.integer)) 
    rownames(occ) <- rownames(tmp) 
    
    min_sample <- min(rowSums(occ))
    
    ## Run Rarefaction
    out <- rarefy(occ, min_sample)
    
    sr_obs <- specnumber(occ)
    
    ## Store results
    dt_res <- data.frame(out, sr_obs)
    names(dt_res) <- c("SR_rarefy", "SR_observed")
    
    dt_res$site <- sites[i]
    
    estimators[[i]] <- dt_res
  } 
  
  res <- do.call(rbind, estimators)
  return(res)
}
# Example
#res <- demonRarefy(samples = coms)
