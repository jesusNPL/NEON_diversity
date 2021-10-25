library(tidyverse)
library(ape)
library(picante)

setwd("Dropbox/Macrosystems_NEON/")

comNEON <- readRDS("DATA/Plants_NEON/vegetation_ALL_NEON_clean.rds")

comNEON <- comNEON %>% 
        filter(Year == 2018) %>% 
        select(decimalLatitude, decimalLongitude, elevation, nlcdClass, 
               domainID, siteID, plotID, SciNames, percentCover) %>% 
        mutate(pID = plotID) %>% 
        separate(pID, sep = "_", into = c("site", "pnum")) %>% 
        mutate(pnum = str_remove(pnum, "^0+")) %>% 
        mutate(plotID_SPEC = paste0("P_", pnum))

saveRDS(comNEON, file = "DATA/Plants_NEON/commNEON_LONG.rds")

comNEON <- comNEON %>% 
        select(domainID, siteID, plotID, plotID_SPEC, pnum, SciNames, percentCover)

sites <- unique(comNEON$siteID)
sitesAnna <- unique(x$Site)

sites <- sort(sites[sites %in% sitesAnna])

comNEON <- comNEON %>% 
        filter(siteID %in% sites)

##### Load and clean data #####

### Phylogeny
phyNEON <- read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S3.nex")

### Traits
trtNEON <- read.csv("DATA/Traits/Traits_BIEN/GapFilling/imputedTRAITS.csv")

trtNEON <- trtNEON %>% 
        select(Species, mean_LNCPLDryMass, mean_LArea, mean_LAreaPLDryMass, 
               mean_LDryMass, mean_LCCPLDryMass, mean_LLifeSpan, mean_LCCPLNC, 
               mean_LNCPLArea, mean_LStomatalPLArea, mean_LFreshMass)

### Community data
commNEON <- read.csv("DATA/Plants_NEON/veg_cover_compiled.csv")

commNEON <- commNEON %>% 
        mutate(pID = plotID) %>% 
        separate(pID, sep = "_", into = c("Site", "pnum")) %>% 
        mutate(pnum = str_remove(pnum, "^0+")) %>% 
        mutate(plotID_SPEC = paste0("P_", pnum)) %>% 
        select(Site, plotID, plotID_SPEC, pnum, comment, todo, 
               tree_cov, herb_cov, plotyear_herb, plotyear_tree, tot_herb, 
               everything())

commNEON <- commNEON %>% 
        filter(todo == "Maybe Out" | todo == "Keep")
        
saveRDS(commNEON, file = "DATA/Plants_NEON/veg_cover_compiled_clean.rds")

##### Match data by community #####

sitesNEON <- sort(unique(comNEON$siteID))

### Match community and phylogenetic data 

matched <- matchNEON_phy(samples = comNEON, phy = phyNEON, 
                         sites = sitesNEON, nSites = length(sitesNEON))

matched_PhyComm <- matched[[1]]

save(matched_PhyComm, file = "DATA/matchPhyloComm/matchPhyloComm.RData")
save(matched, file = "DATA/matchPhyloComm/matches_&_samples.RData")

### Match community and trait data

matchedTrait <- matched_TraitComm(samples = comNEON, traits = trtNEON, 
                                  sites = sitesNEON, nSites = length(sitesNEON))

save(matchedTrait, file = "DATA/matchTraitComm/matchTraitComm.RData")


