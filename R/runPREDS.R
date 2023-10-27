##### PREPARE DATA ######

### Phylogenetic dimension

load("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_q0.RData")
phyQ0 <- list(res$fits[[1]], res$fits[[2]], res$fits[[3]], res$fits[[4]], res$fits[[10]])

load("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_q1.RData")
phyQ1 <- list(res$fits[[10]])

load("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_q2.RData")
phyQ2 <- list(res$fits[[10]])

phyQ <- append(phyQ0, phyQ1)
phyQ <- append(phyQ, phyQ2)

rm(res, phyQ0, phyQ1, phyQ2)

metricPhylo <- c("PD", "PDz", "MPD", "MPDz", "qDTM_p", "qDTM_p", "qDTM_p")
metricSpec <- c("SD", "SD", "MSD", "MSD", "qDSM_s", "qDSM_s", "qDSM_s")

qs <- c("qNA", "qNA", "qNA", "qNA", "q0", "q1", "q2")

save(phyQ, metricPhylo, metricSpec, qs, 
     file = "output/Predictions/data/phy_data.RData")

### Trait dimension
load("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_q0.RData")
trtQ0 <- list(res$fits[[1]], res$fits[[2]], res$fits[[3]], res$fits[[9]])

load("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_q1.RData")
trtQ1 <- list(res$fits[[9]])

load("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_q2.RData")
trtQ2 <- list(res$fits[[9]])

trtQ <- append(trtQ0, trtQ1)
trtQ <- append(trtQ, trtQ2)

metricTrait <- c("MFD", "MFDz", "SR", "qDTM_f", "qDTM_f", "qDTM_f")
metricSpecT <- c("MSD", "MSD", "MSD", "qDSM_s", "qDSM_s", "qDSM_s")

qsT <- c("qNA", "qNA", "qNA", "q0", "q1", "q2")


save(trtQ, metricTrait, metricSpecT, qsT, 
     file = "output/Predictions/data/trt_data.RData")

### Taxonomic dimension
load("Results/Regressions/taxo-spec/scale/reg_Taxo_Spec_Alpha_THRESH.RData")
taxQ <- list(res$fits[[1]], res$fits[[3]], res$fits[[5]])

metricTaxo <- c("SRich_g", "Shannon_g", "Simpson_g")
metricSpectTaxo <- c("S_s", "H_s", "D_g")

qsTax <- c("q0", "q1", "q2")

save(taxQ, metricTaxo, metricSpectTaxo, qsTax, 
     file = "output/Predictions/data/tax_data.RData")

### Spectral species phylogeny
load("Results/Regressions/NEW/phylo_SS.RData")

metricPhylo_SS <- c("MPD", "PD", "q0DTM_p", "q1DTM_p", "q2DTM_p")
metricSpec_SS <- c("SRich_s", "SRich_s", "SRich_s", "Shannon_s", "Simpson_s")

qs_SS <- c("q0", "q0", "q0", "q1", "q2")

### Spectral species trait 
load("Results/Regressions/NEW/trait_SS.RData")

metricTrait_SS <- c("MFD", "q0DTM_f", "q1DTM_f", "q2DTM_f")
metricSpecT_SS <- c("SRich_s", "SRich_s", "Shannon_s", "Simpson_s")

qsT_SS <- c("q0", "q0", "q1", "q2")

##### RUN PREDICTIONS #####

source("R/NEON_diversity/R/Functions/demonPRED.R")

### Phylogenetic dimension 
load("output/Predictions/data/phy_data.RData")

## Inits
nMetrics <- length(metricPhylo)

resPhylo <- list()

for(i in 1:nMetrics) { 
  
  print(paste0("Run predictions for metric ", metricPhylo[i], " and q order ", qs[i]))
  
  fit <- phyQ[[i]]
  
  res <- demonPrediction(fit = fit, variable = metricPhylo[i], covariable = metricSpec[i], 
                         dimension = "Phylogeny", Q = qs[i], nCores = 20, nIters = 2000)
  
  resPhylo[[i]] <- res
  
}

save(resPhylo, file = "output/Predictions/results/predict_phylo.Rdata")

rm(list = ls())

gc()

### Trait dimension
source("R/NEON_diversity/R/Functions/demonPRED.R")

load("output/Predictions/data/trt_data.RData")

## Inits
metricTrait <- c("MFD", "MFDz", "SR", "qDTM_f", "qDTM_f", "qDTM_f")

nMetrics <- length(metricTrait)

resTrait <- list()

for(i in 1:nMetrics) { 
  
  print(paste0("Run predictions for metric ", metricTrait[i], " and q order ", qsT[i]))
  
  fit <- trtQ[[i]]
  
  res <- demonPrediction(fit = fit, variable = metricTrait[i], covariable = metricSpecT[i], 
                         dimension = "Trait", Q = qsT[i], nCores = 28, nIters = 2000)
  
  resTrait[[i]] <- res
  
}

save(resTrait, file = "output/Predictions/results/predict_trait.Rdata")

rm(list = ls())

gc()

### Taxonomic dimension
source("R/NEON_diversity/R/Functions/demonPRED.R")

load("output/Predictions/data/tax_data.RData")

## Inits
nMetrics <- length(metricTaxo)

resTaxo <- list()

for(i in 1:nMetrics) { 
  
  print(paste0("Run predictions for metric ", metricTaxo[i], " and q order ", qsTax[i]))
  
  fit <- taxQ[[i]]
  
  res <- demonPrediction(fit = fit, variable = metricTaxo[i], covariable = metricSpectTaxo[i], 
                         dimension = "Taxonomy", Q = qsTax[i], nCores = 28, nIters = 2000)
  
  resTaxo[[i]] <- res
  
}

save(resTaxo, file = "output/Predictions/results/predict_taxo.Rdata")

rm(list = ls())

gc()



##### Run Spectral species #####
source("R/NEON_diversity/R/Functions/demonPRED.R")

### Spectral species phylogeny
load("Results/Regressions/NEW/phylo_SS.RData")

metricPhylo_SS <- c("MPD", "PD", "q0DTM_p", "q1DTM_p", "q2DTM_p")
metricSpec_SS <- c("SRich_s", "SRich_s", "SRich_s", "Shannon_s", "Simpson_s")

qs_SS <- c("q0", "q0", "q0", "q1", "q2")

## Inits
nMetrics <- length(metricPhylo_SS)

resPhylo_SS <- list()

for(i in 1:nMetrics) { 
  
  print(paste0("Run predictions for metric ", metricPhylo_SS[i], " and q order ", qs_SS[i]))
  
  fit <- fit_phylo_SS[[i]]
  
  res <- demonPrediction(fit = fit, variable = metricPhylo_SS[i], covariable = metricSpec_SS[i], 
                         dimension = "Phylo_SS", Q = qs_SS[i], nCores = 28, nIters = 2000)
  
  resPhylo_SS[[i]] <- res
  
}

save(resPhylo_SS, file = "output/Predictions/results/predict_phylo_SS.Rdata")

rm(list = ls())

gc()

### Spectral species trait 
source("R/NEON_diversity/R/Functions/demonPRED.R")

load("Results/Regressions/NEW/trait_SS.RData")

metricTrait_SS <- c("MFD", "q0DTM_f", "q1DTM_f", "q2DTM_f")
metricSpecT_SS <- c("SRich_s", "SRich_s", "Shannon_s", "Simpson_s")

qsT_SS <- c("q0", "q0", "q1", "q2")

## Inits
nMetrics <- length(metricTrait_SS)

resTrait_SS <- list()

for(i in 1:nMetrics) { 
  
  print(paste0("Run predictions for metric ", metricTrait_SS[i], " and q order ", qsT_SS[i]))
  
  fit <- fit_traitSS[[i]]
  
  res <- demonPrediction(fit = fit, variable = metricTrait_SS[i], covariable = metricSpecT_SS[i], 
                         dimension = "Trait_SS", Q = qsT_SS[i], nCores = 28, nIters = 2000)
  
  resTrait_SS[[i]] <- res
  
}

save(resTrait_SS, file = "output/Predictions/results/predict_trait_SS.Rdata")

rm(list = ls())

gc()



library(mcp)
library(tidybayes)
library(tidyverse)

options(mc.cores = 24)  # Speed up sampling

load("Results/RegDATA/SAM/phylo_spec_DIS_SAM_NEON.RData")

dat <- phylo_spec_NEON_table_SAM_q0 %>% 
  select(Site, plotID, MPD, MPDz, MSD, qDTM_p, qDTM_s) %>% 
  drop_na()


##### Inits #####
iters = 25000
chains = 4

## Model PD

mod_MPD <- list(
  MPDz ~ 1 + MSD,  # intercept + slope
  1 + (1|Site) ~ 0 + MSD,  # joined slope
  1 ~ 0,  # joined plateau
  1 ~ 1  # disjoined plateau
) 

## Run model
fit_mpd <- mcp(data = dat, 
              model = mod_MPD, 
              chains = chains, 
              iter = iters) 

summary(fit_mpd)

p_MPD <- plot(fit_mpd)
p_MPD

pp_MPD <- predict(fit_mpd)

comb <- cbind(dat[, c(1:3)], pp_MPD)
