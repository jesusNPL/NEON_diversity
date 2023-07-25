##### Run Bayesian regressions for NEON phylo-spec diversity #####
library(tidymodels)

dir.create("Results/Regressions/phylo-spec/splitData")

load("Results/RegDATA/SAM/phylo_spec_DIS_SAM_NEON.RData")

source("R/NEON_diversity/R/Functions/GOD_Bayes_phylo_spec.R")

### Q0
# Split data 70%-30% into training set and test set
data_split_q0 <- phylo_spec_NEON_table_SAM_q0 %>% 
  initial_split(prop = 0.70, strata = "Site")

# Extract data in each split
data_train_q0 <- training(data_split_q0)
data_test_q0 <- testing(data_split_q0)

### Q1
# Split data 70%-30% into training set and test set
data_split_q1 <- phylo_spec_NEON_table_SAM_q1 %>% 
  initial_split(prop = 0.70, strata = "Site")

# Extract data in each split
data_train_q1 <- training(data_split_q1)
data_test_q1 <- testing(data_split_q1)

### Q2
# Split data 70%-30% into training set and test set
data_split_q2 <- phylo_spec_NEON_table_SAM_q2 %>% 
  initial_split(prop = 0.70, strata = "Site")

# Extract data in each split
data_train_q2 <- training(data_split_q2)
data_test_q2 <- testing(data_split_q2)

save(data_split_q0, data_split_q1, data_split_q2, 
     file = "Results/RegDATA/SAM/phylo_spec_DIS_SAM_NEON_split.RData")

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

### Run Bayesian models
#for(j in 1:nQ) { 
  ## Q 0
  reg_q0 <- god_BayReg_phylo(resMetrics = data_train_q0,#phylo_spec_NEON_table_q0, 
                             Q = QS[1], nMetrics = 10, 
                             pathSave = paste0("Results/Regressions/phylo-spec/splitData/reg_Phylo_Spec_DIS_split_", QS[1] , ".RData"), 
                             nChains = 4, 
                             nIters = 5000, 
                             nCores = 4, 
                             control = list(adapt_delta = 0.99),
                             engine = "cmdstanr") 
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


########## ----------------------- FINAL ANALYSES ----------------- ############
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
                 nChains = 4, nIters = 5000, nCores = 24, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q1, 
                 Q = QS[2], nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/reg_PhyloSpec_DIS_SAM_", 
                                   QS[2], ".RData"), scale = FALSE,
                 nChains = 4, nIters = 5000, nCores = 24, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q2, 
                 Q = QS[3], nMetrics = 10,  
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/reg_PhyloSpec_DIS_SAM_", 
                                   QS[3], ".RData"),  scale = FALSE,
                 nChains = 4, nIters = 5000, nCores = 24, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q3, 
                 Q = QS[4], nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/reg_PhyloSpec_DIS_SAM_", 
                                   QS[4], ".RData"),  scale = FALSE,
                 nChains = 4, nIters = 5000, nCores = 24, 
                 control = list(adapt_delta = 0.99), engine = "cmdstanr")

##### regressions with scale #####
god_BayReg_phylo(resMetrics = data_train_q0, #phylo_spec_NEON_table_SAM_q0, 
                 Q = QS[1], 
                 nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_SPLIT_", 
                                   QS[1], ".RData"), 
                 scale = TRUE,
                 nChains = 4, 
                 nIters = 5000, 
                 nCores = 8, 
                 control = list(adapt_delta = 0.99), 
                 engine = "cmdstanr")

god_BayReg_phylo(resMetrics = data_train_q1, #phylo_spec_NEON_table_SAM_q1, 
                 Q = QS[2], 
                 nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_SPLIT_", 
                                   QS[2], ".RData"),  
                 scale = TRUE,
                 nChains = 4, 
                 nIters = 5000, 
                 nCores = 8, 
                 control = list(adapt_delta = 0.99),
                 engine = "cmdstanr")

god_BayReg_phylo(resMetrics = data_train_q2, #phylo_spec_NEON_table_SAM_q2, 
                 Q = QS[3], 
                 nMetrics = 10, 
                 pathSave = paste0("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_SPLIT_", 
                                   QS[3], ".RData"), 
                 scale = TRUE,
                 nChains = 4, 
                 nIters = 5000, 
                 nCores = 8, 
                 control = list(adapt_delta = 0.99), 
                 engine = "cmdstanr")

#god_BayReg_phylo(resMetrics = phylo_spec_NEON_table_SAM_q3, 
 #                Q = QS[4], nMetrics = 10, 
  #               pathSave = paste0("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_", 
   #                                QS[4], ".RData"),  scale = TRUE,
    #             nChains = 4, nIters = 5000, nCores = 24, 
     #            control = list(adapt_delta = 0.99), engine = "cmdstanr")

##### Run regressions by habitat #####

library(tidyverse)
library(brms)

source("R/NEON_diversity/R/Functions/GOD_Bayes_phylo_spec.R")

# data on metrics by NEON plots
load("Results/RegDATA/SAM/phylo_spec_DIS_SAM_NEON.RData")

# NEON distributional data
NEON_hab <- read.csv("Results/RegDATA/NEON_metadata_MATCH_plotID.csv")

phylo_spec_NEON_table_q0 <- left_join(phylo_spec_NEON_table_SAM_q0, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

phylo_spec_NEON_table_q1 <- left_join(phylo_spec_NEON_table_SAM_q1, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

phylo_spec_NEON_table_q2 <- left_join(phylo_spec_NEON_table_SAM_q2, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

phylo_spec_NEON_table_q3 <- left_join(phylo_spec_NEON_table_SAM_q3, NEON_hab, 
                                      by = c("Site" = "site", "plotID"))

### Inits
QS <- c("q0", "q1", "q2", "q3")

nQ <- length(QS) 

habitat <- sort(unique(phylo_spec_NEON_table_q0$nlcdClass))
habitat <- c(habitat[2], habitat[5], habitat[6], habitat[7], habitat[10])
nhabitat <- length(habitat)

dir.create("Results/Regressions/phylo-spec/Habitat/Scale")

### Run Bayesian models
for(j in 1:nhabitat) { 
  print(habitat[j])
  
  ## Q 0
  data_q0 <- phylo_spec_NEON_table_q0 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_phylo(resMetrics = data_q0, 
                   Q = QS[1], nMetrics = 10,
                   pathSave = paste0("Results/Regressions/phylo-spec/Habitat/Scale/reg_Phylo_Spec_DIS_", 
                                     QS[1], "_", habitat[j], ".RData"), scale = TRUE,
                   nChains = 4, nIters = 3000, nCores = 24, 
                   control = list(adapt_delta = 0.99), engine = "cmdstanr") 
  ## Q 1 
  data_q1 <- phylo_spec_NEON_table_q1 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_phylo(resMetrics = data_q1, 
                   Q = QS[2], nMetrics = 10, 
                   pathSave = paste0("Results/Regressions/phylo-spec/Habitat/Scale/reg_Phylo_Spec_DIS_", 
                                     QS[2], "_", habitat[j], ".RData"), scale = TRUE,
                   nChains = 4, nIters = 3000, nCores = 24, 
                   control = list(adapt_delta = 0.99), engine = "cmdstanr") 
  ## Q 2
  data_q2 <- phylo_spec_NEON_table_q2 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_phylo(resMetrics = data_q2, 
                   Q = QS[3], nMetrics = 10, 
                   pathSave = paste0("Results/Regressions/phylo-spec/Habitat/Scale/reg_Phylo_Spec_DIS_", 
                                     QS[3], "_", habitat[j] ,".RData"), scale = TRUE,
                   nChains = 4, nIters = 3000, nCores = 24, 
                   control = list(adapt_delta = 0.99), engine = "cmdstanr") 
  ## Q 3 
  data_q3 <- phylo_spec_NEON_table_q3 %>% 
    filter(nlcdClass == habitat[[j]])
  
  god_BayReg_phylo(resMetrics = data_q3, 
                   Q = QS[4], nMetrics = 10, 
                   pathSave = paste0("Results/Regressions/phylo-spec/Habitat/Scale/reg_Phylo_Spec_DIS_", 
                                     QS[4], "_", habitat[j],".RData"), scale = TRUE,
                   nChains = 4, nIters = 3000, nCores = 24, 
                   control = list(adapt_delta = 0.99), engine = "cmdstanr") 
  
}
