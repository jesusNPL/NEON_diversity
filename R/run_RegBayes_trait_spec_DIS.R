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
                             Q = QS[1], nMetrics = 9, 
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

###################### ------------------- Multilelvel 2 ----------------------------- ######################
  
##### Run Bayesian regressions for NEON trait-spec classic diversity metrics #####
  
load("Results/RegDATA/trait_spec_MET_NEON.RData")
  
source("R/NEON_diversity/R/Functions/GOD_Bayes_trait_spec.R")
  
### Inits
QS <- c("q0", "q1", "q2", "q3")

# Q 0
reg_q0 <- god_BayReg_trait_met_ML2(resMetrics = trait_spec_NEON_table_MET_q0, 
                                   Q = QS[1], nMetrics = 5, 
                                   pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_MET_", QS[1] ,"_ML2.RData"), 
                                   nChains = 4, nIters = 2000, nCores = 24, engine = "cmdstanr")

# Q 1
reg_q1 <- god_BayReg_trait_met_ML2(resMetrics = trait_spec_NEON_table_MET_q1, 
                                   Q = QS[2], nMetrics = 5, 
                                   pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_MET_", QS[2] ,"_ML2.RData"), 
                                   nChains = 4, nIters = 2000, nCores = 24, engine = "cmdstanr")

# Q 2
reg_q2 <- god_BayReg_trait_met_ML2(resMetrics = trait_spec_NEON_table_MET_q2, 
                                   Q = QS[3], nMetrics = 5, 
                                   pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_MET_", QS[3] ,"_ML2.RData"), 
                                   nChains = 4, nIters = 2000, nCores = 24, engine = "cmdstanr")

# Q 3
reg_q3 <- god_BayReg_trait_met_ML2(resMetrics = trait_spec_NEON_table_MET_q3, 
                                   Q = QS[4], nMetrics = 5, 
                                   pathSave = paste0("Results/Regressions/trait-spec/reg_Trait_Spec_MET_", QS[4] ,"_ML2.RData"), 
                                   nChains = 4, nIters = 2000, nCores = 24, engine = "cmdstanr")