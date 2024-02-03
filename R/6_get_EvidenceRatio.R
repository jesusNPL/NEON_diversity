library(brms)
library(tidyverse)

########## ---------- Global models Distance based metrics ---------- ########## 

##### Load fitted models for the trait dimension and run demonER #####
source("R/NEON_diversity/R/Functions/demonER.R")
## Q0
load("Results/Regressions/Reanalyses/trait_DIS/BLM_trait_DIS_alpha_q0_scale.RData")
resq0 <- res$fits
## Q1
load("Results/Regressions/Reanalyses/trait_DIS/BLM_trait_DIS_alpha_q1_scale.RData")
resq1 <- res$fits
## Q2
load("Results/Regressions/Reanalyses/trait_DIS/BLM_trait_DIS_alpha_q2_scale.RData")
resq2 <- res$fits

### Run using alpha 0.05
r_mtd_005 <- demonER(fit = resq0[[1]], metric = "MTD", alpha = 0.05, quantile = "full")
r_mtdz_005 <- demonER(fit = resq0[[2]], metric = "MTDz", alpha = 0.05, quantile = "full")
r_qdtm0_005 <- demonER(fit = resq0[[4]], metric = "qDTM_t0", alpha = 0.05, quantile = "full")
r_qdtm1_005 <- demonER(fit = resq1[[4]], metric = "qDTM_t1", alpha = 0.05, quantile = "full")
r_qdtm2_005 <- demonER(fit = resq2[[4]], metric = "qDTM_t2", alpha = 0.05, quantile = "full")

## Run using alpha 0.1
r_mtd_01 <- demonER(fit = resq0[[1]], metric = "MTD", alpha = 0.1, quantile = "full")
r_mtdz_01 <- demonER(fit = resq0[[2]], metric = "MTDz", alpha = 0.1, quantile = "full")
r_qdtm0_01 <- demonER(fit = resq0[[4]], metric = "qDTM_t0", alpha = 0.1, quantile = "full")
r_qdtm1_01 <- demonER(fit = resq1[[4]], metric = "qDTM_t1", alpha = 0.1, quantile = "full")
r_qdtm2_01 <- demonER(fit = resq2[[4]], metric = "qDTM_t2", alpha = 0.1, quantile = "full")

er_trait_all <- rbind(r_mtd_005, r_mtdz_005, r_qdtm0_005, r_qdtm1_005, r_qdtm2_005, 
                      r_mtd_01, r_mtdz_01, r_qdtm0_01, r_qdtm1_01, r_qdtm2_01) 

er_trait_all <- er_trait_all %>% 
  mutate(Dimension = "Trait") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_trait_all, 
          file = "output/EviRatio/Reanalyses/ER_trait_Distance.csv")

##### Load fitted models for the phylogenetic dimension and run demonER #####
source("R/NEON_diversity/R/Functions/demonER.R")

## Q0
load("Results/Regressions/Reanalyses/phylo_DIS/BLM_phylo_DIS_alpha_q0_scale.RData")
resq0 <- res$fits
## Q1
load("Results/Regressions/Reanalyses/phylo_DIS/BLM_phylo_DIS_alpha_q1_scale.RData")
resq1 <- res$fits
## Q2
load("Results/Regressions/Reanalyses/phylo_DIS/BLM_phylo_DIS_alpha_q2_scale.RData")
resq2 <- res$fits

### Run using alpha 0.05
r_pd_005 <- demonER(fit = resq0[[1]], metric = "PD", alpha = 0.05, quantile = "full")
r_pdz_005 <- demonER(fit = resq0[[2]], metric = "PDz", alpha = 0.05, quantile = "full")
r_mpd_005 <- demonER(fit = resq0[[3]], metric = "MPD", alpha = 0.05, quantile = "full")
r_mpdz_005 <- demonER(fit = resq0[[4]], metric = "MPDz", alpha = 0.05, quantile = "full")
r_qdtm0_005 <- demonER(fit = resq0[[10]], metric = "qDPM_p0", alpha = 0.05, quantile = "full")
r_qdtm1_005 <- demonER(fit = resq1[[10]], metric = "qDPM_p1", alpha = 0.05, quantile = "full")
r_qdtm2_005 <- demonER(fit = resq2[[10]], metric = "qDPM_p2", alpha = 0.05, quantile = "full")

## Run using alpha 0.1
r_pd_01 <- demonER(fit = resq0[[1]], metric = "PD", alpha = 0.1, quantile = "full")
r_pdz_01 <- demonER(fit = resq0[[2]], metric = "PDz", alpha = 0.1, quantile = "full")
r_mpd_01 <- demonER(fit = resq0[[3]], metric = "MPD", alpha = 0.1, quantile = "full")
r_mpdz_01 <- demonER(fit = resq0[[4]], metric = "MPDz", alpha = 0.1, quantile = "full")
r_qdtm0_01 <- demonER(fit = resq0[[10]], metric = "qDPM_p0", alpha = 0.1, quantile = "full")
r_qdtm1_01 <- demonER(fit = resq1[[10]], metric = "qDPM_p1", alpha = 0.1, quantile = "full")
r_qdtm2_01 <- demonER(fit = resq2[[10]], metric = "qDPM_p2", alpha = 0.1, quantile = "full")

er_phylo_all <- rbind(r_pd_005, r_pdz_005, r_mpd_005, r_mpdz_005, r_qdtm0_005, r_qdtm1_005, r_qdtm2_005, 
                      r_pd_01, r_pdz_01, r_mpd_01, r_mpdz_01, r_qdtm0_01, r_qdtm1_01, r_qdtm2_01)

er_phylo_all <- er_phylo_all %>% 
  mutate(Dimension = "Phylogeny") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_phylo_all, 
          file = "output/EviRatio/Reanalyses/ER_phylogeny_Distance.csv")

### Combine both sets
er_distance <- bind_rows(er_trait_all, 
                         er_phylo_all)

write_csv(er_distance, 
          file = "output/EviRatio/Reanalyses/ER_trait_phylogeny_DISTANCE.csv")

########## ---------- Quantile models Distance based metrics ---------- ########## 

##### Quantile models for the trait dimension and run demonER #####

### Run ER for quantile models
source("R/NEON_diversity/R/Functions/demonER.R")

quantiles <- c(0.05, 0.1, 0.25, 0.50, 0.75, 0.90, 0.95)
nQuant <- length(quantiles)
metrics <- c("MTD", "MTDz", "qDTM_t0", "qDTM_t1", "qDTM_t2")

### Load data 

## MTD
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_MTD_0595.RData")

## MTDz
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_MTDz_0595.RData")

# qDTM 0
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_qDTM0_0595.RData")

## qDTM 1
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_qDTM1_0595.RData")

## qDTM 2
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_qDTM2_0595.RData")

### Run ER for quantiles 
r_mtd <- list()
r_mtdz <- list()
r_qdtm0 <- list()
r_qdtm1<- list()
r_qdtm2 <- list()

r_mtd_01 <- list()
r_mtdz_01 <- list()
r_qdtm0_01 <- list()
r_qdtm1_01 <- list()
r_qdtm2_01 <- list()

for(i in 1:length(quantiles)) {
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.05"))
  
  r_mtd[[i]] <- demonER(fit = quantREG_mtd[[i]], metric = "MTD", alpha = 0.05, 
                        quantile = quantiles[i]) 
  r_mtdz[[i]] <- demonER(fit = quantREG_mtdz[[i]], metric = "MTDz", alpha = 0.05, 
                         quantile = quantiles[i])
  r_qdtm0[[i]] <- demonER(fit = quantREG_qdtm0[[i]], metric = "qDTM_t0", alpha = 0.05, 
                          quantile = quantiles[i])
  r_qdtm1[[i]] <- demonER(fit = quantREG_qdtm1[[i]], metric = "qDTM_t1", alpha = 0.05, 
                          quantile = quantiles[i])
  r_qdtm2[[i]] <- demonER(fit = quantREG_qdtm2[[i]], metric = "qDTM_t2", alpha = 0.05, 
                          quantile = quantiles[i])
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.1"))
  
  r_mtd_01[[i]] <- demonER(fit = quantREG_mtd[[i]], metric = "MTD", alpha = 0.1, 
                           quantile = quantiles[i]) 
  r_mtdz_01[[i]] <- demonER(fit = quantREG_mtdz[[i]], metric = "MTDz", alpha = 0.1, 
                            quantile = quantiles[i])  
  r_qdtm0_01[[i]] <- demonER(fit = quantREG_qdtm0[[i]], metric = "qDTM_t0", alpha = 0.1, 
                             quantile = quantiles[i])
  r_qdtm1_01[[i]] <- demonER(fit = quantREG_qdtm1[[i]], metric = "qDTM_t1", alpha = 0.1, 
                             quantile = quantiles[i])
  r_qdtm2_01[[i]] <- demonER(fit = quantREG_qdtm2[[i]], metric = "qDTM_t2", alpha = 0.1, 
                             quantile = quantiles[i])

}

### Combine outputs
r_mtd <- do.call(rbind, r_mtd)
r_mtdz <- do.call(rbind, r_mtdz)
r_qdtm0 <- do.call(rbind, r_qdtm0)
r_qdtm1 <- do.call(rbind, r_qdtm1)
r_qdtm2 <- do.call(rbind, r_qdtm2)

r_mtd_01 <- do.call(rbind, r_mtd_01)
r_mtdz_01 <- do.call(rbind, r_mtdz_01)
r_qdtm0_01 <- do.call(rbind, r_qdtm0_01)
r_qdtm1_01 <- do.call(rbind, r_qdtm1_01)
r_qdtm2_01 <- do.call(rbind, r_qdtm2_01)

er_trait_quantile_distance <- rbind(r_mtd, r_mtdz, r_qdtm0, r_qdtm1, r_qdtm2, 
                                    r_mtd_01, r_mtdz_01, r_qdtm0_01, r_qdtm1_01, r_qdtm2_01)

### Save outputs
er_trait_quantile_distance <- er_trait_quantile_distance %>% 
  mutate(Dimension = "Trait") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_trait_quantile_distance, 
          file = "output/EviRatio/Reanalyses/ER_trait_quantile_Distance.csv")

##### Quantile models for the phylogenetic dimension and run demonER #####

source("R/NEON_diversity/R/Functions/demonER.R")

quantiles <- c(0.05, 0.1, 0.25, 0.50, 0.75, 0.90, 0.95)
nQuant <- length(quantiles)
metrics <- c("PD", "PDz", "MPD", "MPDz", "qDPM_p0", "qDPM_p1", "qDPM_p2")

### Load data 
## PD
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_PD_0595.RData")

## PDz
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_PDz_0595.RData")

## MPD
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_MPD_0595.RData")

## MPDz
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_MPDz_0595.RData")

## qDPM 0
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_qDPM0_0595.RData")

## qDPM 1
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_qDPM1_0595.RData")

## qDPM 2
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_qDPM2_0595.RData")

### Run ER for quantiles 
r_pd <- list()
r_pdz <- list()
r_mpd <- list()
r_mpdz <- list()
r_qdpm0 <- list()
r_qdpm1<- list()
r_qdpm2 <- list()

r_pd_01 <- list() 
r_pdz_01 <- list()
r_mpd_01 <- list() 
r_mpdz_01 <- list()
r_qdpm0_01 <- list()
r_qdpm1_01 <- list()
r_qdpm2_01 <- list()

for(i in 1:length(quantiles)) {
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.05"))
  
  r_pd[[i]] <- demonER(fit = quantREG_pd[[i]], metric = "PD", alpha = 0.05, 
                       quantile = quantiles[i])
  r_pdz[[i]] <- demonER(fit = quantREG_pdz[[i]], metric = "PDz", alpha = 0.05, 
                        quantile = quantiles[i])
  r_mpd[[i]] <- demonER(fit = quantREG_mpd[[i]], metric = "MPD", alpha = 0.05, 
                        quantile = quantiles[i])
  r_mpdz[[i]] <- demonER(fit = quantREG_mpdz[[i]], metric = "MPDz", alpha = 0.05, 
                         quantile = quantiles[i]) 
  r_qdpm0[[i]] <- demonER(fit = quantREG_qdpm0[[i]], metric = "qDPM_p0", alpha = 0.05, 
                          quantile = quantiles[i])
  r_qdpm1[[i]] <- demonER(fit = quantREG_qdpm1[[i]], metric = "qDPM_p1", alpha = 0.05, 
                          quantile = quantiles[i])
  r_qdpm2[[i]] <- demonER(fit = quantREG_qdpm2[[i]], metric = "qDPM_p2", alpha = 0.05, 
                          quantile = quantiles[i])
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.1"))
  
  r_pd_01[[i]] <- demonER(fit = quantREG_pd[[i]], metric = "PD", alpha = 0.1, 
                          quantile = quantiles[i])
  r_pdz_01[[i]] <- demonER(fit = quantREG_pdz[[i]], metric = "PDz", alpha = 0.1, 
                           quantile = quantiles[i])
  r_mpd_01[[i]] <- demonER(fit = quantREG_mpd[[i]], metric = "MPD", alpha = 0.1, 
                           quantile = quantiles[i])
  r_mpdz_01[[i]] <- demonER(fit = quantREG_mpdz[[i]], metric = "MPDz", alpha = 0.1, 
                            quantile = quantiles[i]) 
  r_qdpm0_01[[i]] <- demonER(fit = quantREG_qdpm0[[i]], metric = "qDPM_p0", alpha = 0.1, 
                             quantile = quantiles[i])
  r_qdpm1_01[[i]] <- demonER(fit = quantREG_qdpm1[[i]], metric = "qDPM_p1", alpha = 0.1, 
                             quantile = quantiles[i])
  r_qdpm2_01[[i]] <- demonER(fit = quantREG_qdpm2[[i]], metric = "qDPM_p2", alpha = 0.1, 
                             quantile = quantiles[i])
  
}

### Combine outputs
r_pd <- do.call(rbind, r_pd)
r_pdz <- do.call(rbind, r_pdz)
r_mpd <- do.call(rbind, r_mpd)
r_mpdz <- do.call(rbind, r_mpdz)
r_qdpm0 <- do.call(rbind, r_qdpm0)
r_qdpm1 <- do.call(rbind, r_qdpm1)
r_qdpm2 <- do.call(rbind, r_qdpm2)

r_pd_01 <- do.call(rbind, r_pd_01)
r_pdz_01 <- do.call(rbind, r_pdz_01)
r_mpd_01 <- do.call(rbind, r_mpd_01)
r_mpdz_01 <- do.call(rbind, r_mpdz_01)
r_qdpm0_01 <- do.call(rbind, r_qdpm0_01)
r_qdpm1_01 <- do.call(rbind, r_qdpm1_01)
r_qdpm2_01 <- do.call(rbind, r_qdpm2_01)

er_phylo_quantile_distance <- rbind(r_pd, r_pdz, r_mpd, r_mpdz, r_qdpm0, r_qdpm1, r_qdpm2, 
                                    r_pd_01, r_pdz_01, r_mpd_01, r_mpdz_01, r_qdpm0_01, r_qdpm1_01, r_qdpm2_01)

### Save outputs
er_phylo_quantile_distance <- er_phylo_quantile_distance %>% 
  mutate(Dimension = "Phylogeny") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_phylo_quantile_distance, 
          file = "output/EviRatio/Reanalyses/ER_phylo_quantile_Distance.csv")

########## ---------- Global models Spectral Species metrics ---------- ########## 
source("R/NEON_diversity/R/Functions/demonER.R")

##### Taxonomic dimension #####

### Load models
load("Results/Regressions/Reanalyses/taxo_SS/BLM_taxo_SS_alpha.RData")

### Run using alpha 0.05
r_s_005 <- demonER(fit = res$fits[[1]], metric = "S_q0", alpha = 0.05, quantile = "full")
r_h_005 <- demonER(fit = res$fits[[2]], metric = "H_q1", alpha = 0.05, quantile = "full")
r_d_005 <- demonER(fit = res$fits[[3]], metric = "D_q2", alpha = 0.05, quantile = "full")

## Run using alpha 0.1
r_s_01 <- demonER(fit = res$fits[[1]], metric = "S_q0", alpha = 0.1, quantile = "full")
r_h_01 <- demonER(fit = res$fits[[2]], metric = "H_q1", alpha = 0.1, quantile = "full")
r_d_01 <- demonER(fit = res$fits[[3]], metric = "D_q2", alpha = 0.1, quantile = "full")

### Combine and save outputs
er_taxonomy_SS <- rbind(r_s_005, r_h_005, r_d_005, 
                        r_s_01, r_h_01, r_d_01)

er_taxonomy_SS <- er_taxonomy_SS %>% 
  mutate(Dimension = "Taxonomy") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_taxonomy_SS, 
          file = "output/EviRatio/Reanalyses/ER_taxonomy_SS.csv")

##### Trait dimension #####

### Load models 
load("Results/Regressions/Reanalyses/taxo_SS/BLM_trait_SS_alpha.RData")

### Run using alpha 0.05
r_qdtm0_005 <- demonER(fit = res$fits[[1]], metric = "qDTM_p0", alpha = 0.05, quantile = "full")
r_qdtm1_005 <- demonER(fit = res$fits[[2]], metric = "qDTM_p1", alpha = 0.05, quantile = "full")
r_qdtm2_005 <- demonER(fit = res$fits[[3]], metric = "qDTM_p2", alpha = 0.05, quantile = "full")

## Run using alpha 0.1
r_qdtm0_01 <- demonER(fit = res$fits[[1]], metric = "qDTM_p0", alpha = 0.1, quantile = "full")
r_qdtm1_01 <- demonER(fit = res$fits[[2]], metric = "qDTM_p1", alpha = 0.1, quantile = "full")
r_qdtm2_01 <- demonER(fit = res$fits[[3]], metric = "qDTM_p2", alpha = 0.1, quantile = "full")

### Combine and save results
er_trait_SS <- rbind(r_qdtm0_005, r_qdtm1_005, r_qdtm2_005, 
                     r_qdtm0_01, r_qdtm1_01, r_qdtm2_01)

er_trait_SS <- er_trait_SS %>% 
  mutate(Dimension = "Trait") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_trait_SS, 
          file = "output/EviRatio/Reanalyses/ER_trait_SS.csv")

##### Phylogenetic dimension #####

### Load models 
load("Results/Regressions/Reanalyses/taxo_SS/BLM_phylo_SS_alpha.RData")

### Run using alpha 0.05
r_qdpm0_005 <- demonER(fit = res$fits[[1]], metric = "qDPM_p0", alpha = 0.05, quantile = "full")
r_qdpm1_005 <- demonER(fit = res$fits[[2]], metric = "qDPM_p1", alpha = 0.05, quantile = "full")
r_qdpm2_005 <- demonER(fit = res$fits[[3]], metric = "qDPM_p2", alpha = 0.05, quantile = "full")

## Run using alpha 0.1
r_qdpm0_01 <- demonER(fit = res$fits[[1]], metric = "qDPM_p0", alpha = 0.1, quantile = "full")
r_qdpm1_01 <- demonER(fit = res$fits[[2]], metric = "qDPM_p1", alpha = 0.1, quantile = "full")
r_qdpm2_01 <- demonER(fit = res$fits[[3]], metric = "qDPM_p2", alpha = 0.1, quantile = "full")

### Combine and save results
er_phylogeny_SS <- rbind(r_qdpm0_005, r_qdpm1_005, r_qdpm2_005, 
                         r_qdpm0_01, r_qdpm1_01, r_qdpm2_01)

er_phylogeny_SS <- er_phylogeny_SS %>% 
  mutate(Dimension = "Phylogeny") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_phylogeny_SS, 
          file = "output/EviRatio/Reanalyses/ER_phylogeny_SS.csv")

##### Combine all dimensions #####
er_SS <- bind_rows(er_taxonomy_SS, 
                   er_trait_SS, 
                   er_phylogeny_SS)

write_csv(er_SS, 
          file = "output/EviRatio/Reanalyses/ER_taxonomy_trait_phylogeny_SS.CSV")


########## ---------- Quantile models Spectral Species metrics ---------- ########## 

##### Taxonomic dimension #####

### Inits
source("R/NEON_diversity/R/Functions/demonER.R")
quantiles <- c(0.05, 0.1, 0.25, 0.50, 0.75, 0.90, 0.95)
nQuant <- length(quantiles)
metrics <- c("S_q0", "H_q1", "D_q2")

### SR ~ SSR
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/quant_taxo_SR_SS_0595.RData")

### H ~ SSH
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/quant_taxo_H_SS_0595.RData")

### D ~ SSD
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/quant_taxo_D_SS_0595.RData")

### Run ER for quantiles 
r_sr <- list()
r_h <- list()
r_d <- list()

for(i in 1:nQuant) {
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.05"))
  
  r_sr[[i]] <- demonER(fit = quantREG_sr[[i]], metric = "S_q0", alpha = 0.05, 
                       quantile = quantiles[i])
  r_h[[i]] <- demonER(fit = quantREG_h[[i]], metric = "H_q1", alpha = 0.05, 
                      quantile = quantiles[i])
  r_d[[i]] <- demonER(fit = quantREG_d[[i]], metric = "D_q2", alpha = 0.05, 
                      quantile = quantiles[i])
  
}

### Combine outputs
r_sr <- do.call(rbind, r_sr)
r_h <- do.call(rbind, r_h)
r_d <- do.call(rbind, r_d)

### Save outputs
er_taxo_quantile <- rbind(r_sr, r_h, r_d)

er_taxo_quantile <- er_taxo_quantile %>% 
  mutate(Dimension = "Taxonomy") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_taxo_quantile, 
          file = "output/EviRatio/Reanalyses/ER_taxonomy_quantile_SS.csv")

##### Phylogenetic dimension #####  

### Inits
source("R/NEON_diversity/R/Functions/demonER.R")
quantiles <- c(0.05, 0.1, 0.25, 0.50, 0.75, 0.90, 0.95)
nQuant <- length(quantiles)
metrics <- c("qDPM_q0", "qDPM_q1", "qDPM_q2")

### qDPM_q0
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/phylo_q0DPM_SS_quantile.RData")

### qDPM_q1
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/phylo_q1DPM_SS_quantile.RData")

### qDPM_q2
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/phylo_q2DPM_SS_quantile.RData")

### get ER for quantiles
q0dpm_ss <- list()
q1dpm_ss <- list()
q2dpm_ss <- list()

for(i in 1:nQuant) { 
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.05"))
  
  q0dpm_ss[[i]] <- demonER(fit = fit_q0DPM_SS[[i]], metric = metrics[1], alpha = 0.05, 
                           quantile = quantiles[i]) 
  q1dpm_ss[[i]] <- demonER(fit = fit_q1DPM_SS[[i]], metric = metrics[2], alpha = 0.05, 
                           quantile = quantiles[i])
  q2dpm_ss[[i]] <- demonER(fit = fit_q2DPM_SS[[i]], metric = metrics[3], alpha = 0.05, 
                           quantile = quantiles[i])
  
}

### Combine outputs
q0dpm_ss <- do.call(rbind, q0dpm_ss)
q1dpm_ss <- do.call(rbind, q1dpm_ss)
q2dpm_ss <- do.call(rbind, q2dpm_ss)

er_phylo_SS_quantile <- rbind(q0dpm_ss, q1dpm_ss, q2dpm_ss)

### Save outputs
er_phylo_SS_quantile <- er_phylo_SS_quantile %>% 
  mutate(Dimension = "Phylogeny") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_phylo_SS_quantile, 
          file = "output/EviRatio/Reanalyses/ER_phylogeny_quantile_SS.csv")

##### Trait dimension #####  

### Inits
source("R/NEON_diversity/R/Functions/demonER.R")
quantiles <- c(0.05, 0.1, 0.25, 0.50, 0.75, 0.90, 0.95)
nQuant <- length(quantiles)
metrics <- c("qDTM_q0", "qDTM_q1", "qDTM_q2")

### qDTM_q0
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/trait_q0DTM_SS_quantile.RData")

### qDTM_q1
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/trait_q1DTM_SS_quantile.RData")

### qDTM_q2
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/trait_q2DTM_SS_quantile.RData")

### get ER for quantiles
q0dtm_ss <- list()
q1dtm_ss <- list()
q2dtm_ss <- list()

for(i in 1:nQuant) { 
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.05"))
  
  q0dtm_ss[[i]] <- demonER(fit = fit_q0DTM_SS[[i]], metric = metrics[1], alpha = 0.05, 
                           quantile = quantiles[i]) 
  q1dtm_ss[[i]] <- demonER(fit = fit_q1DTM_SS[[i]], metric = metrics[2], alpha = 0.05, 
                           quantile = quantiles[i])
  q2dtm_ss[[i]] <- demonER(fit = fit_q2DTM_SS[[i]], metric = metrics[3], alpha = 0.05, 
                           quantile = quantiles[i])
  
}

### Combine outputs
q0dtm_ss <- do.call(rbind, q0dtm_ss)
q1dtm_ss <- do.call(rbind, q1dtm_ss)
q2dtm_ss <- do.call(rbind, q2dtm_ss)

er_trait_SS_quantile <- rbind(q0dtm_ss, q1dtm_ss, q2dtm_ss)

### Save outputs
er_trait_SS_quantile <- er_trait_SS_quantile %>% 
  mutate(Dimension = "Trait") %>% 
  select(Dimension, Metric, Quantile, Alpha, Parameters, everything())

write_csv(er_trait_SS_quantile, 
          file = "output/EviRatio/Reanalyses/ER_trait_quantile_SS.csv")

##### Combine all dimensions #####

### All dimensions in one table
er_quantile_SS <- bind_rows(er_taxo_quantile, 
                            er_trait_SS_quantile, 
                            er_phylo_SS_quantile)

write_csv(er_quantile_SS, 
          file = "output/EviRatio/Reanalyses/ER_taxonomy_trait_phylogeny_quantile_SS.csv")
