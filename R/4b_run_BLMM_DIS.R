library(tidyverse)
# Auxiliary function 
source("R/NEON_diversity/R/Functions/GOD_Bayes_phylo_spec.R")

########## ---------- Phylogenetic dimension ---------- ########## 

##### Load data #####
## Q0
phylo_spec_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q0.csv")
## Q1
phylo_spec_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q1.csv")
## Q2
phylo_spec_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q2.csv")

### Inits
# Diversity order
QS <- c("q0", "q1", "q2") 
# Number of metrics
nMetrics <- 10
# Number of generations
nIters <- 5000 
# Burnin - set to 20%
burnin <- nIters/5 
# Number of NUTS
nChains <- 4 
# Number of cores
nCores <- 20
# Backend
engine <- "cmdstanr" 
# Path to save the results
saveFilePath <- "Results/Regressions/Reanalyses/phylo_DIS/"

##### Run BLMM between distance metrics #####

### Q0 BLM not scaling spectra metrics
god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q0, 
                 Q = QS[1], 
                 nMetrics = nMetrics, 
                 pathSave = paste0(saveFilePath, "BLM_phylo_DIS_alpha_", QS[1], ".RData"), 
                 scale = FALSE,
                 nChains = nChains, 
                 nIters = nIters, 
                 burnin = burnin, 
                 nCores = nCores, 
                 control = list(adapt_delta = 0.99), 
                 engine = "cmdstanr")

### Q0 BLM scaling spectra metrics
god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q0, 
                 Q = QS[1], 
                 nMetrics = nMetrics, 
                 pathSave = paste0(saveFilePath, "BLM_phylo_DIS_alpha_", QS[1], "_scale.RData"), 
                 scale = TRUE,
                 nChains = nChains, 
                 nIters = nIters, 
                 burnin = burnin, 
                 nCores = nCores, 
                 control = list(adapt_delta = 0.99), 
                 engine = "cmdstanr")

### Q1 BLM not scaling spectra metrics
god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q1, 
                 Q = QS[2], 
                 nMetrics = nMetrics, 
                 pathSave = paste0(saveFilePath, "BLM_phylo_DIS_alpha_", QS[2], ".RData"), 
                 scale = FALSE,
                 nChains = nChains, 
                 nIters = nIters, 
                 burnin = burnin, 
                 nCores = nCores, 
                 control = list(adapt_delta = 0.99), 
                 engine = "cmdstanr")

### Q1 BLM scaling spectra metrics
god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q1, 
                 Q = QS[2], 
                 nMetrics = nMetrics, 
                 pathSave = paste0(saveFilePath, "BLM_phylo_DIS_alpha_", QS[2], "_scale.RData"), 
                 scale = TRUE,
                 nChains = nChains, 
                 nIters = nIters, 
                 burnin = burnin, 
                 nCores = nCores, 
                 control = list(adapt_delta = 0.99), 
                 engine = "cmdstanr")

### Q2 BLM not scaling spectra metrics
god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q2, 
                 Q = QS[3], 
                 nMetrics = nMetrics, 
                 pathSave = paste0(saveFilePath, "BLM_phylo_DIS_alpha_", QS[3], ".RData"), 
                 scale = FALSE,
                 nChains = nChains, 
                 nIters = nIters, 
                 burnin = burnin, 
                 nCores = nCores, 
                 control = list(adapt_delta = 0.99), 
                 engine = "cmdstanr")

### Q2 BLM scaling spectra metrics
god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q2, 
                 Q = QS[3], 
                 nMetrics = nMetrics, 
                 pathSave = paste0(saveFilePath, "BLM_phylo_DIS_alpha_", QS[3], "_scale.RData"), 
                 scale = TRUE,
                 nChains = nChains, 
                 nIters = nIters, 
                 burnin = burnin, 
                 nCores = nCores, 
                 control = list(adapt_delta = 0.99), 
                 engine = "cmdstanr")

########## ---------- Trait dimension ---------- ########## 

# Auxiliary function 
source("R/NEON_diversity/R/Functions/GOD_Bayes_trait_spec.R")

##### Load data #####
## Q0
trait_spec_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q0.csv") %>% 
  mutate(MSD3_spec = MSD_spec) %>% 
  select(siteID, plotID, # Group level
         MFD_trait, MFDz_trait, SR_trait, qDTM_trait, # trait diversity metrics
         MSD_spec, MSD2_spec, MSD3_spec, qDTM_spec) # spectral diversity metrics
  
## Q1
trait_spec_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q1.csv") %>% 
  mutate(MSD3_spec = MSD_spec) %>% 
  select(siteID, plotID, # Group level
         MFD_trait, MFDz_trait, SR_trait, qDTM_trait, # trait diversity metrics
         MSD_spec, MSD2_spec, MSD3_spec, qDTM_spec) # spectral diversity metrics

## Q2
trait_spec_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q2.csv") %>% 
  mutate(MSD3_spec = MSD_spec) %>% 
  select(siteID, plotID, # Group level
         MFD_trait, MFDz_trait, SR_trait, qDTM_trait, # trait diversity metrics
         MSD_spec, MSD2_spec, MSD3_spec, qDTM_spec) # spectral diversity metrics

### Inits
# Diversity order
QS <- c("q0", "q1", "q2") 
# Number of metrics
nMetrics <- 4
# Number of generations
nIters <- 5000 
# Burnin - set to 20%
burnin <- nIters/5 
# Number of NUTS
nChains <- 4 
# Number of cores
nCores <- 20
# Backend
engine <- "cmdstanr" 
# Path to save the results
saveFilePath <- "Results/Regressions/Reanalyses/trait_DIS/"

##### Run BLMM between distance metrics ##### 

## Q0 scaling spectra metrics
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q0, 
                     Q = QS[1], 
                     scaled = TRUE, 
                     pathSave = paste0(saveFilePath, "BLM_trait_DIS_alpha_", QS[1], "_scale.RData"), 
                     nMetrics = nMetrics, 
                     nIters = nIters, 
                     burnin = burnin, 
                     nChains = nChains, 
                     nCores = nCores, 
                     control = list(adapt_delta = 0.99), 
                     engine = engine
)

## Q0 not scaling spectra metrics
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q0, 
                     Q = QS[1], 
                     scaled = FALSE, 
                     pathSave = paste0(saveFilePath, "BLM_trait_DIS_alpha_", QS[1], ".RData"), 
                     nMetrics = nMetrics, 
                     nIters = nIters, 
                     burnin = burnin, 
                     nChains = nChains, 
                     nCores = nCores, 
                     control = list(adapt_delta = 0.99), 
                     engine = engine
) 

## Q1 scaling spectra metrics
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q1, 
                     Q = QS[2], 
                     scaled = TRUE, 
                     pathSave = paste0(saveFilePath, "BLM_trait_DIS_alpha_", QS[2], "_scale.RData"), 
                     nMetrics = nMetrics, 
                     nIters = nIters, 
                     burnin = burnin, 
                     nChains = nChains, 
                     nCores = nCores, 
                     control = list(adapt_delta = 0.99), 
                     engine = engine
) 

## Q1 not scaling spectra metrics
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q1, 
                     Q = QS[2], 
                     scaled = FALSE, 
                     pathSave = paste0(saveFilePath, "BLM_trait_DIS_alpha_", QS[2], ".RData"), 
                     nMetrics = nMetrics, 
                     nIters = nIters, 
                     burnin = burnin, 
                     nChains = nChains, 
                     nCores = nCores, 
                     control = list(adapt_delta = 0.99), 
                     engine = engine
) 

## Q2 scaling spectra metrics
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q2, 
                     Q = QS[3], 
                     scaled = TRUE, 
                     pathSave = paste0(saveFilePath, "BLM_trait_DIS_alpha_", QS[3], "_scale.RData"), 
                     nMetrics = nMetrics, 
                     nIters = nIters, 
                     burnin = burnin, 
                     nChains = nChains, 
                     nCores = nCores, 
                     control = list(adapt_delta = 0.99), 
                     engine = engine
) 

## Q2 not scaling spectra metrics
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q2, 
                     Q = QS[3], 
                     scaled = FALSE, 
                     pathSave = paste0(saveFilePath, "BLM_trait_DIS_alpha_", QS[3], ".RData"), 
                     nMetrics = nMetrics, 
                     nIters = nIters, 
                     burnin = burnin, 
                     nChains = nChains, 
                     nCores = nCores, 
                     control = list(adapt_delta = 0.99), 
                     engine = engine
) 
