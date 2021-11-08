library(tidyverse)

source("R/NEON_diversity/R/Functions/GOD_Bayes_taxo_spec.R")

load("C:/Users/jpintole/Dropbox/Macrosystems_NEON/Results/taxoDiversity/taxo_spec_ALPHA_NEON.RData")

##### Run Bayesian regressions for Alpha taxonomic metrics #####

reg_AlphaOBS_taxo <- god_BayReg_alpha_taxo(resMetrics = taxo_spec_obs_NEON_table, nMetrics = 13, 
                                           nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                                           pathSave = "Results/Regressions/taxo-spec/reg_Taxo_Spec_Alpha_OBS.RData")

reg_AlphaTHRESH_taxo <- god_BayReg_alpha_taxo(resMetrics = taxo_spec_thresh_NEON_table, nMetrics = 13, 
                                              nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                                              pathSave = "Results/Regressions/taxo-spec/reg_Taxo_Spec_Alpha_THRESH.RData")

##### Run Bayesian regressions for Alpha-beta-Gamma diversity at site level #####

### Inits
divLevel <- c("Alpha", "Beta", "Gamma")
qL <- c("q0_Richness", "q1_Shannon", "q2_Simpson")
nLevels <- length(divLevel)

### Run Bayesian models
for(i in 1:nLevels) { 
  
  print(divLevel[i]) 
  
  
  god_BayReg_ABG_taxo(resMetrics = taxo_spec_COM_obs_NEON_table, nMetrics = 2, 
                      diversityLevel = divLevel[i], qlevel = qL[1], 
                      nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                      pathSave = paste0("Results/Regressions/taxo-spec/reg_taxo_", divLevel[i], "_COM_", qL[1], ".RData")) 
  
  god_BayReg_ABG_taxo(resMetrics = taxo_spec_COM_obs_NEON_table, nMetrics = 2, 
                      diversityLevel = divLevel[i], qlevel = qL[2], 
                      nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                      pathSave = paste0("Results/Regressions/taxo-spec/reg_taxo_", divLevel[i], "_COM_", qL[2], ".RData"))
  
  god_BayReg_ABG_taxo(resMetrics = taxo_spec_COM_obs_NEON_table, nMetrics = 2, 
                      diversityLevel = divLevel[i], qlevel = qL[3], 
                      nChains = 4, nIters = 5000, nCores = 26, engine = "cmdstanr", 
                      pathSave = paste0("Results/Regressions/taxo-spec/reg_taxo_", divLevel[i], "_COM", qL[3], ".RData"))
  
  
  }

