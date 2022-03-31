library(tidyverse)

source("R/NEON_diversity/R/Functions/GOD_Bayes_taxo_spec.R")

#load("C:/Users/jpintole/Dropbox/Macrosystems_NEON/Results/taxoDiversity/taxo_spec_ALPHA_NEON.RData")

load("Results/RegDATA/taxo_spec_ALPHA_NEON.RData")
load("Results/RegDATA/taxo_spec_BETA_NEON.RData")
##### Run Bayesian regressions for Alpha taxonomic metrics #####

god_BayReg_alpha_taxo(resMetrics = taxo_spec_obs_NEON_table, nMetrics = 13, scale = TRUE, 
                      nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                      pathSave = "Results/Regressions/taxo-spec/scale/reg_Taxo_Spec_Alpha_OBS.RData")

god_BayReg_alpha_taxo(resMetrics = taxo_spec_thresh_NEON_table, nMetrics = 13, scale = TRUE, 
                      nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                      pathSave = "Results/Regressions/taxo-spec/scale/reg_Taxo_Spec_Alpha_THRESH.RData")

##### Run Bayesian regressions for Alpha-beta-Gamma diversity at site level #####

### Inits
divLevel <- c("Alpha", "Beta", "Gamma")
qL <- c("q0_Richness", "q1_Shannon", "q2_Simpson")
nLevels <- length(divLevel)

### Run Bayesian models
for(i in 1:nLevels) { 
  
  print(divLevel[i]) 
  
  
  god_BayReg_ABG_taxo(resMetrics = taxo_spec_COM_thresh_NEON_table, nMetrics = 2, 
                      diversityLevel = divLevel[i], qlevel = qL[1], scale = FALSE, 
                      nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                      pathSave = paste0("Results/Regressions/taxo-spec/reg_taxo_", divLevel[i], "_COM_thresh_", qL[1], ".RData")) 
  
  god_BayReg_ABG_taxo(resMetrics = taxo_spec_COM_thresh_NEON_table, nMetrics = 2, 
                      diversityLevel = divLevel[i], qlevel = qL[2], scale = FALSE , 
                      nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                      pathSave = paste0("Results/Regressions/taxo-spec/reg_taxo_", divLevel[i], "_COM_thresh_", qL[2], ".RData"))
  
  god_BayReg_ABG_taxo(resMetrics = taxo_spec_COM_thresh_NEON_table, nMetrics = 2, 
                      diversityLevel = divLevel[i], qlevel = qL[3], scale = FALSE, 
                      nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                      pathSave = paste0("Results/Regressions/taxo-spec/reg_taxo_", divLevel[i], "_COM_thresh_", qL[3], ".RData"))
  
  
  }

##### Run regressions by habitat #####

library(tidyverse)
library(brms)

source("R/NEON_diversity/R/Functions/GOD_Bayes_taxo_spec.R")

#load("C:/Users/jpintole/Dropbox/Macrosystems_NEON/Results/taxoDiversity/taxo_spec_ALPHA_NEON.RData")

load("Results/RegDATA/taxo_spec_ALPHA_NEON.RData")

# NEON distributional data
NEON_hab <- read.csv("Results/RegDATA/NEON_metadata_MATCH_plotID.csv")

taxo_spec_obs_NEON_table_habitat <- left_join(taxo_spec_obs_NEON_table, 
                                              NEON_hab, 
                                              by = c("Site" = "site", "plotID"))

taxo_spec_thresh_NEON_table_habitat <- left_join(taxo_spec_thresh_NEON_table, 
                                                 NEON_hab, 
                                                 by = c("Site" = "site", "plotID"))

##### Run Bayesian regressions for Alpha taxonomic metrics ##### 
source("R/NEON_diversity/R/Functions/GOD_Bayes_taxo_spec.R")

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

habitat <- sort(unique(taxo_spec_obs_NEON_table_habitat$nlcdClass))
habitat <- c(habitat[2], habitat[5], habitat[6], habitat[7], habitat[10])
nhabitat <- length(habitat)

dir.create("Results/Regressions/taxo-spec/Habitat")

### Run Bayesian models
for(j in 1:nhabitat) { 
  print(habitat[j])
  
  ## Observed
  data_obs <- taxo_spec_obs_NEON_table_habitat %>% 
    filter(nlcdClass == habitat[[j]]) 
  
  god_BayReg_alpha_taxo(resMetrics = data_obs, nMetrics = 13, scale = FALSE, 
                      nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr", 
                      pathSave = paste0("Results/Regressions/taxo-spec/Habitat/reg_Taxo_Spec_Alpha_OBS_", 
                                        habitat[j], ".RData"))

  data_thresh <- taxo_spec_thresh_NEON_table_habitat %>% 
    filter(nlcdClass == habitat[[j]]) 
  
  ## Thresh
  god_BayReg_alpha_taxo(resMetrics = data_thresh, nMetrics = 13, scale = FALSE, 
                      nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr", 
                      pathSave = paste0("Results/Regressions/taxo-spec/Habitat/reg_Taxo_Spec_Alpha_THRESH_", 
                                        habitat[j], ".RData"))

} 
