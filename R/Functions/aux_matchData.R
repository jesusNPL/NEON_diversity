# Function that perform matching data between communities and phylogeny.
# The function returns a list of matched object with the corresponding NEON site name

### Arguments 
# samples: data.frame with NEON plant data
# phy: phylogenetic data
# sites: vector of NEON site names
# nSites: numeric values of length of sites

matchNEON_phy <- function(samples, phy, sites, nSites) { 
  library(picante)
  
  sampNEON <- list()
  #commNEON <- list()
  matchNEON <- list()
  
  for(i in 1:nSites) { 
    
    samp <- samples %>% 
      filter(siteID == sites[i]) 
    # sample data
    sampNEON[[i]] <- samp
    # community data
    comms <- sample2matrix(samp[, c("plotID_SPEC", "percentCover", "SciNames")])
    #commNEON[[i]] <- comms
    # matched phylo-comm
    matched <- match.phylo.comm(phy = phy, comm = comms)
    matchNEON[[i]] <- matched
    
  } 
  
  names(sampNEON) <- sites
  names(matchNEON) <- sites
  
  res <- list(matcPhyComm = matchNEON, samplesNEON = sampNEON) 
  
  return(res)
  
  print("Match phylogenetic and community data, DONE!")
  
}

# Function that perform match between traits and plant data from NEON

### Arguments 
# samples: data.frame with NEON plant data
# phy: trait data
# sites: vector of NEON site names
# nSites: numeric values of length of sites

matched_TraitComm <- function(samples, traits, sites, nSites) {
  
  library(picante)
  
  matchTrait <- list()
  matchComm <- list()
  
  for(j in 1:nSites) { 
    
    samp <- samples %>% 
      filter(siteID == sites[j]) 
    
    sppSamp <- unique(samp$SciNames)
    # trait data
    traitSamp <- traits[traits$species %in% sppSamp, ]
    spptrait <- traitSamp$species
    
    matchTrait[[j]] <- traitSamp
    
    # community data
    sampNew <- samp[samp$SciNames %in% spptrait, ]
    
    comms <- sample2matrix(sampNew[, c("plotID_SPEC", "percentCover", "SciNames")])
    matchComm[[j]] <- comms
    
  } 
  
  names(matchTrait) <- sites 
  names(matchComm) <- sites
  
  res <- list(traitNEON = matchTrait, commNEON = matchComm) 
  
  return(res) 
  
  print("Match trait and community data, DONE!")
  
}

# this function is for match data using data provided by Anna
matchPhyComm_ANNA <- function(coms, phy, sites, nSites) {
  
  matches <- list() 
  headers <- list()
  
  for(i in 1:nSites) { 
    comtpm <- coms %>% 
      filter(Site == sites[i]) 
    header <- comtpm[, 1:11] 
    comm <- comtpm[, 12:ncol(comtpm)] 
    
    matches[[i]] <- picante::match.phylo.comm(phy = phy, 
                                              comm = data.frame(comm)) 
    headers[[i]] <- header
  }
  
  print(paste0("Match phylo-comm Done!")) 
  
  res <- list(matches, headers)
  return(res)
}

