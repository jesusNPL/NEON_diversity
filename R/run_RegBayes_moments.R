##### Run Bayesian regressions for NEON moments #####

load("Results/RegDATA/moments_ALL_NEON.RData")

source("R/NEON_diversity/R/Functions/GOD_Bayes_moments.R")

##### Make regressions #####
god_BayReg_moments(momentsRES = tbl_moments_NEON, 
                   nChains = 4, 
                   nIters = 10000, 
                   nCores = 3, 
                   engine = "cmdstanr", 
                   nMetrics = 6, 
                   pathSaveTaxo = "Results/Regressions/moments/reg_taxoSpec_moments.RData", 
                   pathSaveTrait = "Results/Regressions/moments/reg_traitSpec_moments.RData", 
                   pathSavePhylo = "Results/Regressions/moments/reg_phyloSpec_moments.RData")

##### Make regressions ML2 #####
god_BayReg_moments_ML2(momentsRES = tbl_moments_NEON, 
                   nChains = 4, 
                   nIters = 10000, 
                   nCores = 3, 
                   engine = "cmdstanr", 
                   nMetrics = 6, 
                   pathSaveTaxo = "Results/Regressions/moments/ML2/reg_taxoSpec_moments_ml2.RData", 
                   pathSaveTrait = "Results/Regressions/moments/ML2/reg_traitSpec_moments_ml2.RData", 
                   pathSavePhylo = "Results/Regressions/moments/ML2/reg_phyloSpec_moments_ml2.RData") 
