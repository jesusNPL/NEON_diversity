##### Run Bayesian regressions for NEON phylo-spec diversity #####

load("Results/RegDATA/phylo_spec_DIS_NEON.RData")

source("R/NEON_diversity/R/Functions/GOD_Bayes_phylo_spec.R")

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

### Run Bayesian models
#for(j in 1:nQ) { 
  ## Q 0
  reg_q0 <- god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q0, 
                             Q = QS[1], nMetrics = 10, 
                             pathSave = paste0("Results/Regressions/phylo-spec/reg_Phylo_Spec_DIS_", QS[1] ,".RData"), 
                             nChains = 4, nIters = 5000, nCores = 4, engine = "cmdstanr") 
  ## Q 1
  reg_q1 <- god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q1, 
                             Q = QS[2], nMetrics = 10, 
                             pathSave = paste0("Results/Regressions/phylo-spec/reg_Phylo_Spec_DIS_", QS[2] ,".RData"), 
                             nChains = 4, nIters = 5000, nCores = 4, engine = "cmdstanr") 
  ## Q 2
  reg_q2 <- god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q2, 
                             Q = QS[3], nMetrics = 10, 
                             pathSave = paste0("Results/Regressions/phylo-spec/reg_Phylo_Spec_DIS_", QS[3] ,".RData"), 
                             nChains = 4, nIters = 5000, nCores = 4, engine = "cmdstanr") 
  ## Q 3
  reg_q3 <- god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_q3, 
                             Q = QS[4], nMetrics = 10, 
                             pathSave = paste0("Results/Regressions/phylo-spec/reg_Phylo_Spec_DIS_", QS[4] ,".RData"), 
                             nChains = 4, nIters = 5000, nCores = 4, engine = "cmdstanr")
  
#}

##### Run Bayesian regression by habitat #####
library(tidyverse)
library(brms)
source("R/NEON_diversity/R/Functions/GOD_Bayes_phylo_spec.R")

load("Results/RegDATA/phylo_spec_DIS_NEON.RData")

NEON_hab <- read.csv("Results/RegDATA/NEON_metadata_MATCH_plotID.csv")

phylo_spec_NEON_table_q0 <- left_join(phylo_spec_NEON_table_q0, NEON_hab, 
                  by = c("Site" = "site", "plotID"))

phylo_spec_NEON_table_q1 <- left_join(phylo_spec_NEON_table_q1, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

phylo_spec_NEON_table_q2 <- left_join(phylo_spec_NEON_table_q2, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

phylo_spec_NEON_table_q3 <- left_join(phylo_spec_NEON_table_q3, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

habitat <- sort(unique(phylo_spec_NEON_table_q0$nlcdClass))
habitat <- c(habitat[2], habitat[5], habitat[6], habitat[7], habitat[10])
nhabitat <- length(habitat)

dir.create("Results/Regressions/phylo-spec/Habitat")

### Run Bayesian models
for(j in 1:nhabitat) { 
  print(habitat[j])
  
## Q 0
  data_q0 <- phylo_spec_NEON_table_q0 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_phylo(resMetrics = data_q0, 
                   Q = QS[1], nMetrics = 10,
                   pathSave = paste0("Results/Regressions/phylo-spec/Habitat/reg_Phylo_Spec_DIS_", QS[1], 
                                     "_", habitat[j], ".RData"), 
                   nChains = 4, nIters = 3000, nCores = 24, 
                   control = list(adapt_delta = 0.99), engine = "cmdstanr") 
## Q 1 
  data_q1 <- phylo_spec_NEON_table_q1 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_phylo(resMetrics = data_q1, 
                   Q = QS[2], nMetrics = 10, 
                   pathSave = paste0("Results/Regressions/phylo-spec/Habitat/reg_Phylo_Spec_DIS_", QS[2], 
                                     "_", habitat[j], ".RData"), 
                   nChains = 4, nIters = 3000, nCores = 24, 
                   control = list(adapt_delta = 0.99), engine = "cmdstanr") 
## Q 2
  data_q2 <- phylo_spec_NEON_table_q2 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_phylo(resMetrics = data_q2, 
                   Q = QS[3], nMetrics = 10, 
                   pathSave = paste0("Results/Regressions/phylo-spec/Habitat/reg_Phylo_Spec_DIS_", QS[3], 
                                     "_", habitat[j] ,".RData"), 
                   nChains = 4, nIters = 3000, nCores = 24, 
                   control = list(adapt_delta = 0.99), engine = "cmdstanr") 
## Q 3 
  data_q3 <- phylo_spec_NEON_table_q3 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_phylo(resMetrics = data_q3, 
                   Q = QS[4], nMetrics = 10, 
                   pathSave = paste0("Results/Regressions/phylo-spec/Habitat/reg_Phylo_Spec_DIS_", QS[4], 
                                     "_", habitat[j],".RData"), 
                   nChains = 4, nIters = 3000, nCores = 24, 
                   control = list(adapt_delta = 0.99), engine = "cmdstanr") 

}


########## RUN Bayesian regressions using SAM #############

load("Results/RegDATA/SAM/phylo_spec_DIS_SAM_NEON.RData")

source("R/NEON_diversity/R/Functions/GOD_Bayes_phylo_spec.R")

### Inits
QS <- c("q0", "q1", "q2", "q3")
nQ <- length(QS)

dir.create("Results/Regressions/phylo-spec/SAM/Scale")

##### Regressions with no scale #####
god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q0, 
                 Q = QS[1], nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/reg_PhyloSpec_DIS_SAM_", 
                 QS[1], ".RData"), scale = FALSE,
                 nChains = 4, nIters = 5000, nCores = 8, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q1, 
                 Q = QS[2], nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/reg_PhyloSpec_DIS_SAM_", 
                                   QS[2], ".RData"), scale = FALSE,
                 nChains = 4, nIters = 5000, nCores = 8, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q2, 
                 Q = QS[3], nMetrics = 10,  
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/reg_PhyloSpec_DIS_SAM_", 
                                   QS[3], ".RData"),  scale = FALSE,
                 nChains = 4, nIters = 5000, nCores = 8, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q3, 
                 Q = QS[4], nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/reg_PhyloSpec_DIS_SAM_", 
                                   QS[4], ".RData"),  scale = FALSE,
                 nChains = 4, nIters = 5000, nCores = 8, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

##### regressions with scale #####
god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q0, 
                 Q = QS[1], nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_", 
                                   QS[1], ".RData"), scale = TRUE,
                 nChains = 4, nIters = 5000, nCores = 8, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q1, 
                 Q = QS[2], nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_", 
                                   QS[2], ".RData"),  scale = TRUE,
                 nChains = 4, nIters = 5000, nCores = 8, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q2, 
                 Q = QS[3], nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_", 
                                   QS[3], ".RData"), scale = TRUE,
                 nChains = 4, nIters = 5000, nCores = 8, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q3, 
                 Q = QS[4], nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_", 
                                   QS[4], ".RData"),  scale = TRUE,
                 nChains = 4, nIters = 5000, nCores = 8, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")
