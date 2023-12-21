### Libraries
library(tidyverse)
library(brms)

options(mc.cores = 80)

########## ---------- Phylogenetic dimension ---------- ########## 

##### Load data #####
## Q0
phylo_spec_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q0.csv")
## Q1
phylo_spec_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q1.csv")
## Q2
phylo_spec_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q2.csv")

### Inits
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
# Quantiles 
quants <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
# Number of quantiles
nQuants <- length(quants)

##### Run BQLMM distance metrics - Phylogenetic dimension #####

##### PD and PDz ##### 

## list of results
quantREG_pd <- list()
quantREG_pdz <- list()

### Run BQLM
for(i in 1:nQuants) { 
  
  print(paste0("Starting quantile regression for ", quants[i], " %...")) 
  
  ## PD
  quantREG_pd[[i]] <- brm(formula = bf(PD_phylo ~ scale(SD_spec) + (1|siteID), 
                                       quantile = quants[i]), 
                          data = phylo_spec_NEON_table_q0, 
                          family = asym_laplace(), 
                          iter = nIters, 
                          warmup = burnin, 
                          cores = nCores, 
                          control = list(adapt_delta = 0.99), 
                          backend = engine)
                          
  ## PDz
  quantREG_pdz[[i]] <- brm(formula = bf(PDz_phylo ~ scale(SD_spec) + (1|siteID), 
                                        quantile = quants[i]), 
                           data = phylo_spec_NEON_table_q0, 
                           family = asym_laplace(), 
                           iter = nIters, 
                           warmup = burnin, 
                           cores = nCores, 
                           control = list(adapt_delta = 0.99), 
                           backend = engine)
  
}

### Save model fits
## PD
save(quantREG_pd, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_PD_0595.RData")

## PDz
save(quantREG_pdz, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_PDz_0595.RData")

##### MPD and MPDz #####

## list of results
quantREG_mpd <- list()
quantREG_mpdz <- list()

### Run BQLM
for(i in 1:nQuants) { 
  
  print(paste0("Starting quantile regression for ", quants[i], " %...")) 
  
  ## MPD
  quantREG_mpd[[i]] <- brm(formula = bf(MPD_phylo ~ scale(MSD_spec) + (1|siteID), 
                                        quantile = quants[i]), 
                           data = phylo_spec_NEON_table_q0, 
                           family = asym_laplace(), 
                           iter = nIters, 
                           warmup = burnin, 
                           cores = nCores, 
                           control = list(adapt_delta = 0.99), 
                           backend = engine) 
  
  ## MPDz
  quantREG_mpdz[[i]] <- brm(formula = bf(MPDz_phylo ~ scale(MSD_spec) + (1|siteID), 
                                         quantile = quants[i]), 
                            data = phylo_spec_NEON_table_q0, 
                            family = asym_laplace(), 
                            iter = nIters, 
                            warmup = burnin, 
                            cores = nCores, 
                            control = list(adapt_delta = 0.99), 
                            backend = engine)
  
}

### Save model fits
## MPD
save(quantREG_mpd, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_MPD_0595.RData")

## MPDz
save(quantREG_mpdz, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_MPDz_0595.RData") 

##### qDPM ##### 

## list of results
quantREG_qdpm0 <- list()
quantREG_qdpm1 <- list()
quantREG_qdpm2 <- list()

for(i in 1:nQuants) { 
  
  print(paste0("Starting quantile regression for ", quants[i], " %...")) 
  
  
  ## Q0
  quantREG_qdpm0[[i]] <- brm(formula = bf(qDTM_phylo ~ scale(qDTM_spec) + (1|siteID), 
                                          quantile = quants[i]), 
                             data = phylo_spec_NEON_table_q0, 
                             family = asym_laplace(), 
                             iter = nIters, 
                             warmup = burnin, 
                             cores = nCores, 
                             control = list(adapt_delta = 0.99), 
                             backend = engine)
  
  ## Q1
  quantREG_qdpm1[[i]] <- brm(formula = bf(qDTM_phylo ~ scale(qDTM_spec) + (1|siteID), 
                                          quantile = quants[i]), 
                             data = phylo_spec_NEON_table_q1, 
                             family = asym_laplace(), 
                             iter = nIters, 
                             warmup = burnin, 
                             cores = nCores, 
                             control = list(adapt_delta = 0.99), 
                             backend = engine)
  
  ## Q2
  quantREG_qdpm2[[i]] <- brm(formula = bf(qDTM_phylo ~ scale(qDTM_spec) + (1|siteID), 
                                          quantile = quants[i]), 
                             data = phylo_spec_NEON_table_q2, 
                             family = asym_laplace(), 
                             iter = nIters, 
                             warmup = burnin, 
                             cores = nCores, 
                             control = list(adapt_delta = 0.99), 
                             backend = engine)
  
}

### Save model fits
## qDPM0
save(quantREG_qdpm0, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_qDPM0_0595.RData")

## qDPM1
save(quantREG_qdpm1, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_qDPM1_0595.RData")

## qDPM2
save(quantREG_qdpm2, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_qDPM2_0595.RData")

rm(list = ls())

gc()

########## ---------- Trait dimension ---------- ##########

##### Load data #####
## Q0
trait_spec_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q0.csv") 
## Q1
trait_spec_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q1.csv") 
## Q2
trait_spec_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q2.csv") 

### Inits
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
# Quantiles 
quants <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
# Number of quantiles
nQuants <- length(quants)

##### Run BQLMM distance metrics - trait dimension ##### 

##### MTD, MTDz and SR #####

## list of results
quantREG_mtd <- list()
quantREG_mtdz <- list()
quantREG_sr <- list()

### Run BQLM
for(i in 1:nQuants) { 
  
  print(paste0("Starting quantile regression for ", quants[i], " %...")) 
  
  ## MTD
  quantREG_mtd[[i]] <- brm(formula = bf(MFD_trait ~ scale(MSD_spec) + (1|siteID), 
                                        quantile = quants[i]), 
                           data = trait_spec_NEON_table_q0, 
                           family = asym_laplace(), 
                           iter = nIters, 
                           warmup = burnin, 
                           cores = nCores, 
                           control = list(adapt_delta = 0.99), 
                           backend = engine) 
  
  ## MTDz
  quantREG_mtdz[[i]] <- brm(formula = bf(MFDz_trait ~ scale(MSD_spec) + (1|siteID), 
                                         quantile = quants[i]), 
                            data = trait_spec_NEON_table_q0, 
                            family = asym_laplace(), 
                            iter = nIters, 
                            warmup = burnin, 
                            cores = nCores, 
                            control = list(adapt_delta = 0.99), 
                            backend = engine) 
  
  ## SR
  quantREG_sr[[i]] <- brm(formula = bf(SR_trait ~ scale(MSD_spec) + (1|siteID), 
                                         quantile = quants[i]), 
                            data = trait_spec_NEON_table_q0, 
                            family = asym_laplace(), 
                            iter = nIters, 
                            warmup = burnin, 
                            cores = nCores, 
                            control = list(adapt_delta = 0.99), 
                            backend = engine)
  
}

### Save model fits
## MTD
save(quantREG_mtd, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_MTD_0595.RData")

## MTDz
save(quantREG_mtdz, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_MTDz_0595.RData") 

## SR
save(quantREG_sr, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_SR_0595.RData") 

##### qDTM ##### 

## list of results
quantREG_qdtm0 <- list()
quantREG_qdtm1 <- list()
quantREG_qdtm2 <- list()

for(i in 1:nQuants) { 
  
  print(paste0("Starting quantile regression for ", quants[i], " %...")) 
  
  
  ## Q0
  quantREG_qdtm0[[i]] <- brm(formula = bf(qDTM_trait ~ scale(qDTM_spec) + (1|siteID), 
                                          quantile = quants[i]), 
                             data = trait_spec_NEON_table_q0, 
                             family = asym_laplace(), 
                             iter = nIters, 
                             warmup = burnin, 
                             cores = nCores, 
                             control = list(adapt_delta = 0.99), 
                             backend = engine)
  
  ## Q1
  quantREG_qdtm1[[i]] <- brm(formula = bf(qDTM_trait ~ scale(qDTM_spec) + (1|siteID), 
                                          quantile = quants[i]), 
                             data = trait_spec_NEON_table_q1, 
                             family = asym_laplace(), 
                             iter = nIters, 
                             warmup = burnin, 
                             cores = nCores, 
                             control = list(adapt_delta = 0.99), 
                             backend = engine)
  
  ## Q2
  quantREG_qdtm2[[i]] <- brm(formula = bf(qDTM_trait ~ scale(qDTM_spec) + (1|siteID), 
                                          quantile = quants[i]), 
                             data = trait_spec_NEON_table_q2, 
                             family = asym_laplace(), 
                             iter = nIters, 
                             warmup = burnin, 
                             cores = nCores, 
                             control = list(adapt_delta = 0.99), 
                             backend = engine)
  
}

### Save model fits
## qDTM0
save(quantREG_qdtm0, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_qDTM0_0595.RData")

## qDTM1
save(quantREG_qdtm1, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_qDTM1_0595.RData")

## qDTM2
save(quantREG_qdtm2, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_qDTM2_0595.RData")

rm(list = ls())

gc()

##### End Bayesian quantile models for distance matrices #####