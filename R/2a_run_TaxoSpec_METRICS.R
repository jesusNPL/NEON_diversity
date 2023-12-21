###### ---------- Part 1 - run taxonomic dimension ---------- #####

### Libraries
library(tidyverse)
library(picante)
# remotes::install_github("cran/vegetarian")
library(vegetarian)

### Auxiliary functions
source("R/NEON_diversity/R/Functions/SS_diversity.R")

#### Prepare data #####
# Load plot species composition
load("DATA/matchPhyloComm/matchPhyloComm.RData")

# Site names
sites <- names(matched_PhyComm)

# Number of sites
nSites <- length(sites)

##### Alpha and beta diversity using the ground data #####

## Store results

# Alpha taxonomic diversity at plot level
alphaGround_NEON <- list()

# Alpha, Beta, and Gamma diversity at site level
betaGround_NEON <- list()

### Start calculations 
for (j in 1:nSites) {
  
  print(sites[j])
  
  ### Inits 
  # site community composition
  commGround <- matched_PhyComm[[j]][[2]] 
  # plots within the site
  plotNames <- rownames(commGround)
  # number of plots within the site
  nPlots <- length(plotNames)
  
  ### Run taxonomic diversity
  resSite <- demon_SSppDIV(
    comm = commGround, # plot community composition data
    q1 = 0, # diversity order 0
    q2 = 1, # diversity order 1
    q3 = 2, # diversity order 2
    plotNames = plotNames, # plot names
    nPlots = nPlots, # number of plot in the site 
    site = sites[j] # site name
  )
  
  # Save results at plot level
  alphaGround_NEON[[j]] <- resSite[[1]] 
  # Save results at site level
  betaGround_NEON[[j]] <- resSite[[2]]
  
  print("Alpha-beta diversity for ground data, done!!!")
  
}

##### Save Results #####

### Combine results into data.frames

# Alpha diversity at plot level
alphaGround_NEON_table <- do.call(rbind, alphaGround_NEON) 

# Alpha, Beta, and Gamma at site level
betaGround_NEON_table <- do.call(rbind, betaGround_NEON)

### Save results 
# Plot level results 
alphaGround_NEON_table <- alphaGround_NEON_table %>% 
  select(Site, plotID, everything())

write_csv(alphaGround_NEON_table, 
          file = "Results/taxoDiversity/Reanalyses/taxo_alphaGround_NEON.csv")

save(alphaGround_NEON_table, alphaGround_NEON,
     file = "Results/taxoDiversity/Reanalyses/alphaGround_diversity_NEON.RData"
)

# Site level results 
betaGround_NEON_table <- betaGround_NEON_table %>% 
  select(Site, everything())

write_csv(betaGround_NEON_table, 
          file = "Results/taxoDiversity/Reanalyses/taxo_SiteGround_NEON.csv")

save(betaGround_NEON_table, betaGround_NEON,
     file = "Results/taxoDiversity/Reanalyses/betaGround_diversity_NEON.RData"
)

##### ---------- Part 2 -  Run metrics of taxonomic dimension using CDM based on spectral species ---------- #####

### Libraries
library(tidyverse)
library(vegetarian)

### Auxiliary functions
source("R/NEON_diversity/R/Functions/SS_diversity.R")

#### Prepare data #####
# Load plot species composition
load("DATA/matchPhyloComm/matchPhyloComm.RData")

# Site names
siteNames <- names(matched_PhyComm)

# Number of sites
nSites <- length(siteNames)

### Run diversity metrics using spectral species
# Alpha taxonomic diversity at plot level
alphaSS_NEON <- list()

# Alpha, Beta, and Gamma diversity at site level
betaSS_NEON <- list()

for (j in 1:nSites) {

  # Inits
  print(siteNames[j])

  ### Load CDM based on spectral species
  commSpecThresh <- read_csv(paste0("Results/taxoDiversity/SpectralCommunities/Reanalyses/spectral_CDM_", 
                                    siteNames[j], "_10metrics.csv")) %>% 
    as.data.frame() 
  
  # Assign plotID as rownames
  rownames(commSpecThresh) <- commSpecThresh$plotID 
  
  # Remove plotID
  commSpecThresh <- commSpecThresh %>% 
    select(!plotID)

  # Plot names in the site
  plotNames <- rownames(commSpecThresh) 
  
  # Number of plots in the site
  nPlots <- length(plotNames)

  ### Run metrics of biodiversity using spectral species
  resSite <- demon_SSppDIV(
    comm = commSpecThresh, # CDM based on spectral species 
    q1 = 0, # diversity order 0
    q2 = 1, # diversity order 1 
    q3 = 2, # diversity order 2
    plotNames = plotNames, # plot names
    nPlots = nPlots, # number of plots 
    site = siteNames[j] # siteID
  ) 
  
  # Save results at plot level
  alphaSS_NEON[[j]] <- resSite[[1]] 
  
  # Save results at site level
  betaSS_NEON[[j]] <- resSite[[2]]

  print(paste0("Alpha-beta diversity using SS for site ", siteNames[j],  " done!!!"))

}

### Combine results
alphaSS_NEON_table <- do.call(rbind, alphaSS_NEON)

betaSS_NEON_table <- do.call(rbind, betaSS_NEON)

### Save results
# Plot level results 
alphaSS_NEON_table <- alphaSS_NEON_table %>% 
  select(Site, plotID, everything())

write_csv(alphaSS_NEON_table , 
          file = "Results/taxoDiversity/Reanalyses/taxo_alphaSS_NEON.csv")

save(alphaSS_NEON_table, alphaSS_NEON,
     file = "Results/taxoDiversity/Reanalyses/alphaSS_diversity_NEON.RData"
)

# Site level results 
betaSS_NEON_table <- betaSS_NEON_table %>% 
  select(Site, everything())

write_csv(betaSS_NEON_table, 
          file = "Results/taxoDiversity/Reanalyses/taxo_SiteSS_NEON.csv")

save(betaSS_NEON_table, betaSS_NEON,
     file = "Results/taxoDiversity/Reanalyses/betaSS_diversity_NEON.RData"
)
