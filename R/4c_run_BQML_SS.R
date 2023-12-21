### Libraries
library(tidyverse)
library(brms)

########## ---------- Taxonomic dimension ---------- ##########

##### Load data #####
taxo_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/taxo_spec_SS_NEON.csv") %>% 
  mutate(siteID = Site)

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

##### Run BQLMM spectral species metrics - taxonomic dimension ##### 

##### SR, H, D #####
## list of results 
quantREG_sr <- list()
quantREG_h <- list()
quantREG_d <- list()

for(i in 1:nQuants) { 
  
  print(paste0("Starting quantile regression for ", quants[i], " %...")) 
  
  ## SR
  quantREG_sr[[i]] <- brm(formula = bf(SRich_ground ~ SRich_spec + (1|siteID), 
                                       quantile = quants[i]), 
                          data = taxo_spec_SS_NEON, 
                          family = asym_laplace(), 
                          iter = nIters, 
                          warmup = burnin, 
                          cores = nCores, 
                          control = list(adapt_delta = 0.99), 
                          backend = engine)
  ## Shannon
  quantREG_h[[i]] <- brm(formula = bf(Shannon_ground ~ Shannon_spec + (1|siteID), 
                                      quantile = quants[i]), 
                         data = taxo_spec_SS_NEON, 
                         family = asym_laplace(), 
                         iter = nIters, 
                         warmup = burnin, 
                         cores = nCores, 
                         control = list(adapt_delta = 0.99), 
                         backend = engine)
  ## Simpson
  quantREG_d[[i]] <- brm(formula = bf(Simpson_ground ~ Simpson_spec + (1|siteID), 
                                      quantile = quants[i]), 
                         data = taxo_spec_SS_NEON, 
                         family = asym_laplace(), 
                         iter = nIters, 
                         warmup = burnin, 
                         cores = nCores, 
                         control = list(adapt_delta = 0.99), 
                         backend = engine)
  
}

### Save model fits
## SR
save(quantREG_sr, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/SS/quant_taxo_SR_SS_0595.RData")

## H
save(quantREG_h, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/SS/quant_taxo_H_SS_0595.RData")

## D
save(quantREG_d, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/SS/quant_taxo_D_SS_0595.RData")

rm(list = ls())

gc()

########## ---------- Phylogenetic dimension ---------- ##########

##### Load data #####
phylo_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_SS_NEON.csv") 

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

##### Run BQLMM spectral species metrics - phylogenetic dimension ##### 

##### qDPM ##### 
### List of results
fit_q0DPM_SS <- list()
fit_q1DPM_SS <- list()
fit_q2DPM_SS <- list()

for(i in 1:nQuants) { 
  
  print(paste0("Run quantile regression for theta ", quants[i], " ... patience..."))
  
  ### q0DPM ~ SS richness
  fit_q0DPM_SS[[i]] <- brms::brm(formula = bf(qDTM_phylo_q0 ~ SRich_spec + (1|siteID), 
                                              quantile = quants[i]), 
                                 data = phylo_spec_SS_NEON, 
                                 family = asym_laplace(), 
                                 iter = nIters, 
                                 warmup = burnin, 
                                 cores = nCores, 
                                 control = list(adapt_delta = 0.99), 
                                 backend = engine)
  
  ### q1DPM ~ SS Shannon
  fit_q1DPM_SS[[i]] <- brms::brm(formula = bf(qDTM_phylo_q1 ~ Shannon_spec + (1|siteID), 
                                              quantile = quants[i]), 
                                 data = phylo_spec_SS_NEON, 
                                 family = asym_laplace(), 
                                 iter = nIters, 
                                 warmup = burnin, 
                                 cores = nCores, 
                                 control = list(adapt_delta = 0.99), 
                                 backend = engine)
  
  ### q2DPM ~ SS Simpson
  fit_q2DPM_SS[[i]] <- brms::brm(formula = bf(qDTM_phylo_q2 ~ Simpson_spec + (1|siteID), 
                                              quantile = quants[i]), 
                                 data = phylo_spec_SS_NEON, 
                                 family = asym_laplace(), 
                                 iter = nIters, 
                                 warmup = burnin, 
                                 cores = nCores, 
                                 control = list(adapt_delta = 0.99), 
                                 backend = engine)
  
}

### Save model fits
## q0DPM
save(fit_q0DPM_SS, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/SS/phylo_q0DPM_SS_quantile.RData")

## q1DPM
save(fit_q1DPM_SS, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/SS/phylo_q1DPM_SS_quantile.RData")

## q2DPM
save(fit_q2DPM_SS, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/SS/phylo_q2DPM_SS_quantile.RData")

rm(list = ls())

gc()

########## ---------- Trait dimension ---------- ##########

##### Load data #####
trait_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/trait_spec_SS_NEON.csv") 

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

##### Run BQLMM spectral species metrics - trait dimension ##### 

##### qDTM ##### 
### List of results
fit_q0DTM_SS <- list()
fit_q1DTM_SS <- list()
fit_q2DTM_SS <- list()

for(i in 1:nQuants) { 
  
  print(paste0("Run quantile regression for theta ", quants[i], " ... patience..."))
  
  ### q0DTM ~ SS richness
  fit_q0DTM_SS[[i]] <- brms::brm(formula = bf(qDTM_trait_q0 ~ SRich_spec + (1|siteID), 
                                              quantile = quants[i]), 
                                 data = trait_spec_SS_NEON, 
                                 family = asym_laplace(), 
                                 iter = nIters, 
                                 warmup = burnin, 
                                 cores = nCores, 
                                 control = list(adapt_delta = 0.99), 
                                 backend = engine)
  
  ### q1DTM ~ SS Shannon
  fit_q1DTM_SS[[i]] <- brms::brm(formula = bf(qDTM_trait_q1 ~ Shannon_spec + (1|siteID), 
                                              quantile = quants[i]), 
                                 data = trait_spec_SS_NEON, 
                                 family = asym_laplace(), 
                                 iter = nIters, 
                                 warmup = burnin, 
                                 cores = nCores, 
                                 control = list(adapt_delta = 0.99), 
                                 backend = engine)
  
  ### q2DTM ~ SS Simpson
  fit_q2DTM_SS[[i]] <- brms::brm(formula = bf(qDTM_trait_q2 ~ Simpson_spec + (1|siteID), 
                                              quantile = quants[i]), 
                                 data = trait_spec_SS_NEON, 
                                 family = asym_laplace(), 
                                 iter = nIters, 
                                 warmup = burnin, 
                                 cores = nCores, 
                                 control = list(adapt_delta = 0.99), 
                                 backend = engine)
  
}

### Save model fits
## q0DTM
save(fit_q0DTM_SS, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/SS/trait_q0DTM_SS_quantile.RData")

## q1DTM
save(fit_q1DTM_SS, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/SS/trait_q1DTM_SS_quantile.RData")

## q2DTM
save(fit_q2DTM_SS, 
     file = "Results/Regressions/Reanalyses/QuantileBLMM/SS/trait_q2DTM_SS_quantile.RData")

rm(list = ls())

gc()

##### End Bayesian quantile models for spectral species #####