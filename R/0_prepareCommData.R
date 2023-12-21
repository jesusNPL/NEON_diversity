library(tidyverse)
library(ape)
library(picante)

##### Load and clean plant community composition #####
comNEON <- readRDS("DATA/Plants_NEON/vegetation_ALL_NEON_clean.rds")

comNEON <- comNEON %>% 
        filter(Year == 2018) %>% 
        select(decimalLatitude, decimalLongitude, elevation, nlcdClass, 
               domainID, siteID, plotID, SciNames, percentCover) %>% 
        mutate(pID = plotID) %>% 
        separate(pID, sep = "_", into = c("site", "pnum")) %>% 
        mutate(pnum = str_remove(pnum, "^0+")) %>% 
        mutate(plotID_SPEC = paste0("P_", pnum))

## Save processes data
#saveRDS(comNEON, file = "DATA/Plants_NEON/commNEON_LONG.rds")

##### Match Community composition with phylogeny and traits #####
comNEON <- readRDS("DATA/Plants_NEON/commNEON_LONG.rds")

comNEON <- comNEON %>% 
        select(domainID, siteID, plotID, plotID_SPEC, pnum, SciNames, percentCover)

sites <- unique(comNEON$siteID)

## Select sites with HSI
x <- data.frame(list.files("DATA/Spectra/"))
names(x) <- "header"
xx <- x %>% 
        separate(header, sep = "_", 
                 into = c("Site", "spec", "grd"))

sitesAnna <- unique(xx$Site)

sites <- sort(sites[sites %in% sitesAnna])

comNEON <- comNEON %>% 
        filter(siteID %in% sites)

##### Load and clean data #####

### Phylogeny
phyNEON <- read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S3.nex")

### Traits
trtNEON <- read.csv("DATA/Traits/GapFilling/imputedTRAITS.csv")

# Select leaf traits
trtNEON <- trtNEON %>% 
        select(species, mean_LNCPLDryMass, mean_LArea, mean_LAreaPLDryMass, 
               mean_LDryMass, mean_LCCPLDryMass, mean_LLifeSpan, mean_LCCPLNC, 
               mean_LNCPLArea, mean_LStomatalPLArea, mean_LFreshMass)


##### Match data by community #####
load("R/NEON_diversity/R/Functions/aux_matchData.R")

sitesNEON <- sort(unique(comNEON$siteID))

### Match community and phylogenetic data 

matched <- matchNEON_phy(samples = comNEON, # NEON plots
                         phy = phyNEON, # phylogeny
                         sites = sitesNEON, # site names
                         nSites = length(sitesNEON) # number of sites
                         ) 

matched_PhyComm <- matched[[1]]

save(matched_PhyComm, file = "DATA/matchPhyloComm/matchPhyloComm.RData")
save(matched, file = "DATA/matchPhyloComm/matches_&_samples.RData")

### Match community and trait data

matchedTrait <- matched_TraitComm(samples = comNEON, # NEON plots
                                  traits = trtNEON, # traits 
                                  sites = sitesNEON, # site names
                                  nSites = length(sitesNEON) # number of sites
                                  )

matchedTrait[[1]]$ABBY
matchedTrait[[2]]$ABBY

save(matchedTrait, file = "DATA/matchTraitComm/matchTraitComm_Check.RData")


