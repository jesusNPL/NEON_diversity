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
