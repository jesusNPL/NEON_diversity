##### Run Bayesian regressions for NEON trait-spec distance based diversity metrics #####

load("Results/RegDATA/trait_spec_DIS_NEON.RData")

source("R/NEON_diversity/R/Functions/GOD_Bayes_trait_spec.R")

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

### Run Bayesian models
#for(j in 1:nQ) { 
  ## Q 0
  reg_q0 <- god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q0, 
                             Q = QS[1], nMetrics = 9, scale = TRUE, 
                             pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_DIS_", QS[1] ,".RData"), 
                             nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr") 
  ## Q 1
  reg_q1 <- god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q1, 
                             Q = QS[2], nMetrics = 9, 
                             pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_DIS_", QS[2] ,".RData"), 
                             nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr") 
  ## Q 2
  reg_q2 <- god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q2, 
                             Q = QS[3], nMetrics = 9, 
                             pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_DIS_", QS[3] ,".RData"), 
                             nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr") 
  ## Q 3
  reg_q3 <- god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_q3, 
                             Q = QS[4], nMetrics = 9, 
                             pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_DIS_", QS[4] ,".RData"), 
                             nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr")
  
#}

rm(list = ls())

##### Run Bayesian regressions for NEON trait-spec classic diversity metrics #####

load("Results/RegDATA/trait_spec_MET_NEON.RData")

source("R/NEON_diversity/R/Functions/GOD_Bayes_trait_spec.R")

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

### Run Bayesian models
#for(j in 1:nQ) { 
  ## Q 0
  reg_q0 <- god_BayReg_trait_met(resMetrics = trait_spec_NEON_table_MET_q0, 
                                 Q = QS[1], nMetrics = 5, 
                                 pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_MET_", QS[1] ,".RData"), 
                                 nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr") 
  ## Q 1
  reg_q1 <- god_BayReg_trait_met(resMetrics = trait_spec_NEON_table_MET_q1, 
                                 Q = QS[2], nMetrics = 5, 
                                 pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_MET_", QS[2] ,".RData"), 
                                 nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr") 
  ## Q 2
  reg_q2 <- god_BayReg_trait_met(resMetrics = trait_spec_NEON_table_MET_q2, 
                                 Q = QS[3], nMetrics = 5, 
                                 pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_MET_", QS[3] ,".RData"), 
                                 nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr") 
  ## Q 3
  reg_q3 <- god_BayReg_trait_met(resMetrics = trait_spec_NEON_table_MET_q3, 
                                 Q = QS[4], nMetrics = 5, 
                                 pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_MET_", QS[4] ,".RData"), 
                                 nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr")
  
#}

###################### ------------------- Multilevel 2 ----------------------------- ######################
  
##### Run Bayesian regressions for NEON trait-spec classic diversity metrics #####
  
load("Results/RegDATA/trait_spec_MET_NEON.RData")
  
source("R/NEON_diversity/R/Functions/GOD_Bayes_trait_spec.R")
  
### Inits
QS <- c("q0", "q1", "q2", "q3")

# Q 0
reg_q0 <- god_BayReg_trait_met_ML2(resMetrics = trait_spec_NEON_table_MET_q0, 
                                   Q = QS[1], nMetrics = 5, 
                                   pathSave = paste0("Results/Regressions/trait-spec/ML2/reg_Trait_Spec_MET_", QS[1] ,"_ML2.RData"), 
                                   nChains = 4, nIters = 2000, nCores = 24, engine = "cmdstanr")

# Q 1
reg_q1 <- god_BayReg_trait_met_ML2(resMetrics = trait_spec_NEON_table_MET_q1, 
                                   Q = QS[2], nMetrics = 5, 
                                   pathSave = paste0("Results/Regressions/trait-spec/ML2/reg_Trait_Spec_MET_", QS[2] ,"_ML2.RData"), 
                                   nChains = 4, nIters = 2000, nCores = 24, engine = "cmdstanr")

# Q 2
reg_q2 <- god_BayReg_trait_met_ML2(resMetrics = trait_spec_NEON_table_MET_q2, 
                                   Q = QS[3], nMetrics = 5, 
                                   pathSave = paste0("Results/Regressions/trait-spec/ML2/reg_Trait_Spec_MET_", QS[3] ,"_ML2.RData"), 
                                   nChains = 4, nIters = 2000, nCores = 24, engine = "cmdstanr")

# Q 3
reg_q3 <- god_BayReg_trait_met_ML2(resMetrics = trait_spec_NEON_table_MET_q3, 
                                   Q = QS[4], nMetrics = 5, 
                                   pathSave = paste0("Results/Regressions/trait-spec/ML2/reg_Trait_Spec_MET_", QS[4] ,"_ML2.RData"), 
                                   nChains = 4, nIters = 2000, nCores = 24, engine = "cmdstanr")

##### ML2 distance based metrics
## Q 0
reg_q0 <- god_BayReg_trait_dis_ML2(resMetrics = trait_spec_NEON_table_q0, 
                               Q = QS[1], nMetrics = 9, 
                               pathSave = paste0("Results/Regressions/trait-spec/ML2/reg_Trait_Spec_DIS_", QS[1], "_ML2.RData"), 
                               nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr") 
## Q 1
reg_q1 <- god_BayReg_trait_dis_ML2(resMetrics = trait_spec_NEON_table_q1, 
                               Q = QS[2], nMetrics = 9, 
                               pathSave = paste0("Results/Regressions/trait-spec/ML2/reg_Trait_Spec_DIS_", QS[2], "_ML2.RData"), 
                               nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr") 
## Q 2
reg_q2 <- god_BayReg_trait_dis_ML2(resMetrics = trait_spec_NEON_table_q2, 
                               Q = QS[3], nMetrics = 9, 
                               pathSave = paste0("Results/Regressions/trait-spec/ML2/reg_Trait_Spec_DIS_", QS[3], "_ML2.RData"), 
                               nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr") 
## Q 3
reg_q3 <- god_BayReg_trait_dis_ML2(resMetrics = trait_spec_NEON_table_q3, 
                               Q = QS[4], nMetrics = 9, 
                               pathSave = paste0("Results/Regressions/trait-spec/ML2/reg_Trait_Spec_DIS_", QS[4], "_ML2.RData"), 
                               nChains = 4, nIters = 10000, nCores = 24, engine = "cmdstanr")


#################### ----------- FINAL ANALYSES ------------ ###################
########## RUN Bayesian regressions using SAM #############

load("Results/RegDATA/SAM/trait_spec_DIS_SAM_NEON.RData")

source("R/NEON_diversity/R/Functions/GOD_Bayes_trait_spec.R")

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 
dir.create("Results/Regressions/trait-spec/SAM/Scale")

### Run Bayesian models
#for(j in 1:nQ) { 
## Q 0
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_SAM_q0, 
                     Q = QS[1], nMetrics = 9, scale = TRUE, 
                     pathSave = paste0("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_", 
                                       QS[1] ,".RData"), 
                     nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr") 
## Q 1
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_SAM_q1, 
                     Q = QS[2], nMetrics = 9, scale = TRUE, 
                     pathSave = paste0("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_", 
                                       QS[2] ,".RData"), 
                     nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr") 
## Q 2
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_SAM_q2, 
                     Q = QS[3], nMetrics = 9, scale = TRUE,  
                     pathSave = paste0("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_", 
                                       QS[3] ,".RData"), 
                     nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr") 
## Q 3
god_BayReg_trait_dis(resMetrics = trait_spec_NEON_table_SAM_q3, 
                     Q = QS[4], nMetrics = 9, scale = TRUE,  
                     pathSave = paste0("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_", 
                                       QS[4] ,".RData"), 
                     nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr")

rm(list = ls())

##### Run Bayesian regressions for NEON trait-spec classic diversity metrics #####

load("Results/RegDATA/SAM/trait_spec_MET_SAM_NEON.RData")

source("R/NEON_diversity/R/Functions/GOD_Bayes_trait_spec.R")

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

### Run Bayesian models
#for(j in 1:nQ) { 
## Q 0
god_BayReg_trait_met(resMetrics = trait_spec_NEON_table_MET_SAM_q0, 
                     Q = QS[1], nMetrics = 5, scale = TRUE, 
                     pathSave = paste0("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_MET_SAM_scaled_", QS[1] ,".RData"), 
                     nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr") 
## Q 1
god_BayReg_trait_met(resMetrics = trait_spec_NEON_table_MET_SAM_q1, 
                     Q = QS[2], nMetrics = 5, scale = TRUE, 
                     pathSave = paste0("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_MET_SAM_scaled_", QS[2] ,".RData"), 
                     nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr") 
## Q 2
god_BayReg_trait_met(resMetrics = trait_spec_NEON_table_MET_SAM_q2, 
                     Q = QS[3], nMetrics = 5, scale = TRUE,  
                     pathSave = paste0("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_MET_SAM_scaled_", QS[3] ,".RData"), 
                     nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr") 
## Q 3
god_BayReg_trait_met(resMetrics = trait_spec_NEON_table_MET_SAM_q3, 
                     Q = QS[4], nMetrics = 5, scale = TRUE,  
                     pathSave = paste0("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_MET_SAM_scaled", QS[4] ,".RData"), 
                     nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr")

##### BY HABITATS DISTANCE BASED METRICS #####

library(tidyverse)
library(brms)
source("R/NEON_diversity/R/Functions/GOD_Bayes_trait_spec.R")

load("Results/RegDATA/SAM/trait_spec_DIS_SAM_NEON.RData")

NEON_hab <- read.csv("Results/RegDATA/NEON_metadata_MATCH_plotID.csv")

trait_spec_NEON_table_q0 <- left_join(trait_spec_NEON_table_SAM_q0, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

trait_spec_NEON_table_q1 <- left_join(trait_spec_NEON_table_SAM_q1, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

trait_spec_NEON_table_q2 <- left_join(trait_spec_NEON_table_SAM_q2, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

trait_spec_NEON_table_q3 <- left_join(trait_spec_NEON_table_SAM_q3, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

habitat <- sort(unique(trait_spec_NEON_table_q0$nlcdClass))
habitat <- c(habitat[2], habitat[5], habitat[6], habitat[7], habitat[10])
nhabitat <- length(habitat)

dir.create("Results/Regressions/trait-spec/Habitat")

### Run Bayesian models for distance based metrics
for(j in 1:nhabitat) { 
  
  print(habitat[j])
  
  ## Q 0
  data_q0 <- trait_spec_NEON_table_q0 %>% 
    filter(nlcdClass == habitat[[j]]) 
  
  god_BayReg_trait_dis(resMetrics = data_q0, 
                       Q = QS[1], nMetrics = 9, scaled = TRUE, 
                       pathSave = paste0("Results/Regressions/trait-spec/Habitat/DIS/reg_Trait_Spec_DIS_SAM_scaled_", 
                                         QS[1], "_", habitat[j], ".RData"), 
                       nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr") 
  ## Q 1
  data_q1 <- trait_spec_NEON_table_q1 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_trait_dis(resMetrics = data_q1, 
                       Q = QS[2], nMetrics = 9, scaled = TRUE, 
                       pathSave = paste0("Results/Regressions/trait-spec/Habitat/DIS/reg_Trait_Spec_DIS_SAM_scaled_", 
                                         QS[2], "_", habitat[j], ".RData"), 
                       nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr") 
  ## Q 2
  data_q2 <- trait_spec_NEON_table_q2 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_trait_dis(resMetrics = data_q2, 
                       Q = QS[3], nMetrics = 9, scaled = TRUE,  
                       pathSave = paste0("Results/Regressions/trait-spec/Habitat/DIS/reg_Trait_Spec_DIS_SAM_scaled_", 
                                         QS[3], "_", habitat[j], ".RData"), 
                       nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr") 
  ## Q 3
  data_q3 <- trait_spec_NEON_table_q3 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_trait_dis(resMetrics = data_q3, 
                       Q = QS[4], nMetrics = 9, scaled = TRUE,  
                       pathSave = paste0("Results/Regressions/trait-spec/Habitat/DIS/reg_Trait_Spec_DIS_SAM_scaled_", 
                                         QS[4], "_", habitat[j], ".RData"), 
                       nChains = 4, nIters = 5000, nCores = 24, engine = "cmdstanr")
  
}

gc()


##### BY HABITATS DISTANCE BASED CLASSIC METRICS #####

library(tidyverse)
library(brms)

source("R/NEON_diversity/R/Functions/GOD_Bayes_trait_spec.R")

load("Results/RegDATA/SAM/trait_spec_MET_SAM_NEON.RData")

NEON_hab <- read.csv("Results/RegDATA/NEON_metadata_MATCH_plotID.csv")

trait_spec_NEON_table_MET_q0 <- left_join(trait_spec_NEON_table_MET_SAM_q0, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

trait_spec_NEON_table_MET_q1 <- left_join(trait_spec_NEON_table_MET_SAM_q1, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

trait_spec_NEON_table_MET_q2 <- left_join(trait_spec_NEON_table_MET_SAM_q2, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

trait_spec_NEON_table_MET_q3 <- left_join(trait_spec_NEON_table_MET_SAM_q3, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

habitat <- sort(unique(trait_spec_NEON_table_MET_q0$nlcdClass))
habitat <- c(habitat[2], habitat[5], habitat[6], habitat[7], habitat[10])
nhabitat <- length(habitat)

dir.create("Results/Regressions/trait-spec/Habitat")

### Run Bayesian models for classic based metrics
for(j in 1:nhabitat) { 
  
  print(habitat[j])
  
  ## Q 0
  data_q0 <- trait_spec_NEON_table_MET_q0 %>% 
    filter(nlcdClass == habitat[[j]]) 
  
  god_BayReg_trait_met(resMetrics = data_q0, 
                       Q = QS[1], nMetrics = 5, scaled = TRUE, 
                       pathSave = paste0("Results/Regressions/trait-spec/Habitat/reg_Trait_Spec_MET_SAM_scaled_", 
                                         QS[1], "_", habitat[j], ".RData"),  
                       nChains = 4, nIters = 5000, nCores = 12, engine = "cmdstanr") 
  
  ## Q 1
  data_q1 <- trait_spec_NEON_table_MET_q1 %>% 
    filter(nlcdClass == habitat[[j]]) 
  
  god_BayReg_trait_met(resMetrics = data_q1, 
                       Q = QS[2], nMetrics = 5, scaled = TRUE, 
                       pathSave = paste0("Results/Regressions/trait-spec/Habitat/reg_Trait_Spec_MET_SAM_scaled_", 
                                         QS[2], "_", habitat[j], ".RData"), 
                       nChains = 4, nIters = 5000, nCores = 12, engine = "cmdstanr") 
  
  ## Q 2
  data_q2 <- trait_spec_NEON_table_MET_q2 %>% 
    filter(nlcdClass == habitat[[j]]) 
  
  god_BayReg_trait_met(resMetrics = data_q2, 
                       Q = QS[3], nMetrics = 5, scaled = TRUE,  
                       pathSave = paste0("Results/Regressions/trait-spec/Habitat/reg_Trait_Spec_MET_SAM_scaled_", 
                                         QS[3], "_", habitat[j], ".RData"), 
                       nChains = 4, nIters = 5000, nCores = 12, engine = "cmdstanr") 
  
  ## Q 3
  data_q3 <- trait_spec_NEON_table_MET_q3 %>% 
    filter(nlcdClass == habitat[[j]]) 
  
  god_BayReg_trait_met(resMetrics = data_q3, 
                       Q = QS[4], nMetrics = 5, scaled = TRUE,  
                       pathSave = paste0("Results/Regressions/trait-spec/Habitat/reg_Trait_Spec_MET_SAM_scaled_", 
                                         QS[4], "_", habitat[j], ".RData"), 
                       nChains = 4, nIters = 5000, nCores = 12, engine = "cmdstanr")
  
}
