library(brms)
library(tidyverse)

##### Load fitted models for the trait dimension and run demonER #####

## Q0
load("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_q0.RData")
resq0 <- res$fits
## Q1
load("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_q1.RData")
resq1 <- res$fits
## Q2
load("Results/Regressions/trait-spec/SAM/Scale/reg_Trait_Spec_DIS_SAM_scaled_q2.RData")
resq2 <- res$fits

### Run using alpha 0.05
r_mtd_005 <- demonER(fit = resq0[[2]], metric = "MTDz", alpha = 0.05, quantile = "full")
r_qdtm0_005 <- demonER(fit = resq0[[9]], metric = "qDTM_t0", alpha = 0.05, quantile = "full")
r_qdtm1_005 <- demonER(fit = resq1[[9]], metric = "qDTM_t1", alpha = 0.05, quantile = "full")
r_qdtm2_005 <- demonER(fit = resq2[[9]], metric = "qDTM_t2", alpha = 0.05, quantile = "full")

## Run using alpha 0.1
r_mtd_01 <- demonER(fit = resq0[[2]], metric = "MTDz", alpha = 0.1, quantile = "full")
r_qdtm0_01 <- demonER(fit = resq0[[9]], metric = "qDTM_t0", alpha = 0.1, quantile = "full")
r_qdtm1_01 <- demonER(fit = resq1[[9]], metric = "qDTM_t1", alpha = 0.1, quantile = "full")
r_qdtm2_01 <- demonER(fit = resq2[[9]], metric = "qDTM_t2", alpha = 0.1, quantile = "full")

er_trait_all <- rbind(r_mtd_005, r_qdtm0_005, r_qdtm1_005, r_qdtm2_005, 
                      r_mtd_01, r_qdtm0_01, r_qdtm1_01, r_qdtm2_01)

write.csv(er_trait_all, file = "output/EviRatio/ER_trait_continental.csv")

### Run ER for quantile models
source("R/NEON_diversity/R/Functions/demonER.R")

quantiles <- c(0.05, 0.1, 0.25, 0.50, 0.75, 0.90, 0.95)
nQuant <- length(quantiles)
metrics <- c("MTDz", "qDTM_t0", "qDTM_t1", "qDTM_t2")

# Load data 
load("Results/Regressions/QuantileREG/Trait/quant_trait_MFDz_site_0595.RData")
mtd_0595 <- quantREG_mfd
load("Results/Regressions/QuantileREG/Trait/quant_trait_MFDz_site.RData")
mtd_2575 <- quantREG_mfd

mtd_quant <- list(mtd_0595[[1]], mtd_0595[[2]], mtd_2575[[1]], mtd_2575[[2]], 
                  mtd_2575[[3]], mtd_0595[[3]], mtd_0595[[4]])

load("Results/Regressions/QuantileREG/Trait/quant_trait_qDTM0_site_0595.RData")
qdtm0_0595 <- quantREG_qdtm0
load("Results/Regressions/QuantileREG/Trait/quant_trait_qDTM0_site.RData")
qdtm0_2575 <- quantREG_qdtm0

qdtm0_quant <- list(qdtm0_0595[[1]], qdtm0_0595[[2]], qdtm0_2575[[1]], qdtm0_2575[[2]], 
                    qdtm0_2575[[3]], qdtm0_0595[[3]], qdtm0_0595[[4]])

load("Results/Regressions/QuantileREG/Trait/quant_trait_qDTM1_site_0595.RData")
qdtm1_0595 <- quantREG_qdtm1
load("Results/Regressions/QuantileREG/Trait/quant_trait_qDTM1_site.RData")
qdtm1_2575 <- quantREG_qdtm1

qdtm1_quant <- list(qdtm1_0595[[1]], qdtm1_0595[[2]], qdtm1_2575[[1]], qdtm1_2575[[2]], 
                    qdtm1_2575[[3]], qdtm1_0595[[3]], qdtm1_0595[[4]])

load("Results/Regressions/QuantileREG/Trait/quant_trait_qDTM2_site_0595.RData")
qdtm2_0595 <- quantREG_qdtm2
load("Results/Regressions/QuantileREG/Trait/quant_trait_qDTM2_site.RData")
qdtm2_2575 <- quantREG_qdtm2

qdtm2_quant <- list(qdtm2_0595[[1]], qdtm2_0595[[2]], qdtm2_2575[[1]], qdtm2_2575[[2]], 
                    qdtm2_2575[[3]], qdtm2_0595[[3]], qdtm2_0595[[4]])

### Run ER for quantiles 
r_mtd <- list()
r_qdtm0 <- list()
r_qdtm1<- list()
r_qdtm2 <- list()

r_mtd_01 <- list()
r_qdtm0_01 <- list()
r_qdtm1_01 <- list()
r_qdtm2_01 <- list()

for(i in 1:length(quantiles)) {
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.05"))
  
  r_mtd[[i]] <- demonER(fit = mtd_quant[[i]], metric = "MTDz", alpha = 0.05, quantile = quantiles[i])
  r_qdtm0[[i]] <- demonER(fit = qdtm0_quant[[i]], metric = "qDTM_t0", alpha = 0.05, quantile = quantiles[i])
  r_qdtm1[[i]] <- demonER(fit = qdtm1_quant[[i]], metric = "qDTM_t1", alpha = 0.05, quantile = quantiles[i])
  r_qdtm2[[i]] <- demonER(fit = qdtm2_quant[[i]], metric = "qDTM_t2", alpha = 0.05, quantile = quantiles[i])
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.1"))
  
  r_mtd_01[[i]] <- demonER(fit = mtd_quant[[i]], metric = "MTDz", alpha = 0.1, quantile = quantiles[i])
  r_qdtm0_01[[i]] <- demonER(fit = qdtm0_quant[[i]], metric = "qDTM_t0", alpha = 0.1, quantile = quantiles[i])
  r_qdtm1_01[[i]] <- demonER(fit = qdtm1_quant[[i]], metric = "qDTM_t1", alpha = 0.1, quantile = quantiles[i])
  r_qdtm2_01[[i]] <- demonER(fit = qdtm2_quant[[i]], metric = "qDTM_t2", alpha = 0.1, quantile = quantiles[i])
  
}

r_mtd <- do.call(rbind, r_mtd)
r_qdtm0 <- do.call(rbind, r_qdtm0)
r_qdtm1 <- do.call(rbind, r_qdtm1)
r_qdtm2 <- do.call(rbind, r_qdtm2)

r_mtd_01 <- do.call(rbind, r_mtd_01)
r_qdtm0_01 <- do.call(rbind, r_qdtm0_01)
r_qdtm1_01 <- do.call(rbind, r_qdtm1_01)
r_qdtm2_01 <- do.call(rbind, r_qdtm2_01)

er_trait_quantile <- rbind(r_mtd, r_qdtm0, r_qdtm1, r_qdtm2, 
                           r_mtd_01, r_qdtm0_01, r_qdtm1_01, r_qdtm2_01)

write.csv(er_trait_quantile, file = "output/EviRatio/ER_trait_continental_quantile.csv")

##### Load fitted models for the phylogenetic dimension and run demonER #####
source("R/NEON_diversity/R/Functions/demonER.R")

## Q0
load("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_q0.RData")
resq0 <- res$fits
## Q1
load("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_q1.RData")
resq1 <- res$fits
## Q2
load("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_q2.RData")
resq2 <- res$fits

### Run using alpha 0.05
r_mtd_005 <- demonER(fit = resq0[[4]], metric = "MPDz", alpha = 0.05, quantile = "full")
r_qdtm0_005 <- demonER(fit = resq0[[10]], metric = "qDPM_p0", alpha = 0.05, quantile = "full")
r_qdtm1_005 <- demonER(fit = resq1[[10]], metric = "qDPM_p1", alpha = 0.05, quantile = "full")
r_qdtm2_005 <- demonER(fit = resq2[[10]], metric = "qDPM_p2", alpha = 0.05, quantile = "full")

## Run using alpha 0.1
r_mtd_01 <- demonER(fit = resq0[[4]], metric = "MPDz", alpha = 0.1, quantile = "full")
r_qdtm0_01 <- demonER(fit = resq0[[10]], metric = "qDPM_p0", alpha = 0.1, quantile = "full")
r_qdtm1_01 <- demonER(fit = resq1[[10]], metric = "qDPM_p1", alpha = 0.1, quantile = "full")
r_qdtm2_01 <- demonER(fit = resq2[[10]], metric = "qDPM_p2", alpha = 0.1, quantile = "full")

er_phylo_all <- rbind(r_mtd_005, r_qdtm0_005, r_qdtm1_005, r_qdtm2_005, 
                      r_mtd_01, r_qdtm0_01, r_qdtm1_01, r_qdtm2_01)

write.csv(er_phylo_all, file = "output/EviRatio/ER_phylo_continental.csv")

### Run ER for quantile models
source("R/NEON_diversity/R/Functions/demonER.R")

quantiles <- c(0.05, 0.1, 0.25, 0.50, 0.75, 0.90, 0.95)
nQuant <- length(quantiles)
metrics <- c("MPDz", "qDPM_p0", "qDPM_p1", "qDPM_p2")

# Load data 
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_MPDz_site_0595.RData")
mpd_0595 <- quantREG_mpd
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_MPDz_site.RData")
mpd_2575 <- quantREG_mpd

mpd_quant <- list(mpd_0595[[1]], mpd_0595[[2]], mpd_2575[[1]], mpd_2575[[2]], 
                  mpd_2575[[3]], mpd_0595[[3]], mpd_0595[[4]])

load("Results/Regressions/QuantileREG/Phylo/quant_phylo_qDTM0_site_0595.RData")
qdpm0_0595 <- quantREG_qdtm0
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_qDTM0_site.RData")
qdpm0_2575 <- quantREG_qdtm0

qdpm0_quant <- list(qdpm0_0595[[1]], qdpm0_0595[[2]], qdpm0_2575[[1]], qdpm0_2575[[2]], 
                    qdpm0_2575[[3]], qdpm0_0595[[3]], qdpm0_0595[[4]])

load("Results/Regressions/QuantileREG/Phylo/quant_phylo_qDTM1_site_0595.RData")
qdpm1_0595 <- quantREG_qdtm1
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_qDTM1_site.RData")
qdpm1_2575 <- quantREG_qdtm1

qdpm1_quant <- list(qdpm1_0595[[1]], qdpm1_0595[[2]], qdpm1_2575[[1]], qdpm1_2575[[2]], 
                    qdpm1_2575[[3]], qdpm1_0595[[3]], qdpm1_0595[[4]])

load("Results/Regressions/QuantileREG/Phylo/quant_phylo_qDTM2_site_0595.RData")
qdpm2_0595 <- quantREG_qdtm2
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_qDTM2_site.RData")
qdpm2_2575 <- quantREG_qdtm2

qdpm2_quant <- list(qdpm2_0595[[1]], qdpm2_0595[[2]], qdpm2_2575[[1]], qdpm2_2575[[2]], 
                    qdpm2_2575[[3]], qdpm2_0595[[3]], qdpm2_0595[[4]])

### Run ER for quantiles 
r_mpd <- list()
r_qdpm0 <- list()
r_qdpm1<- list()
r_qdpm2 <- list()

r_mpd_01 <- list()
r_qdpm0_01 <- list()
r_qdpm1_01 <- list()
r_qdpm2_01 <- list()

for(i in 1:length(quantiles)) {
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.05"))
  
  r_mpd[[i]] <- demonER(fit = mpd_quant[[i]], metric = "MPDz", alpha = 0.05, quantile = quantiles[i])
  r_qdpm0[[i]] <- demonER(fit = qdpm0_quant[[i]], metric = "qDPM_p0", alpha = 0.05, quantile = quantiles[i])
  r_qdpm1[[i]] <- demonER(fit = qdpm1_quant[[i]], metric = "qDPM_p1", alpha = 0.05, quantile = quantiles[i])
  r_qdpm2[[i]] <- demonER(fit = qdpm2_quant[[i]], metric = "qDPM_p2", alpha = 0.05, quantile = quantiles[i])
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.1"))
  
  r_mpd_01[[i]] <- demonER(fit = mpd_quant[[i]], metric = "MPDz", alpha = 0.1, quantile = quantiles[i])
  r_qdpm0_01[[i]] <- demonER(fit = qdpm0_quant[[i]], metric = "qDPM_p0", alpha = 0.1, quantile = quantiles[i])
  r_qdpm1_01[[i]] <- demonER(fit = qdpm1_quant[[i]], metric = "qDPM_p1", alpha = 0.1, quantile = quantiles[i])
  r_qdpm2_01[[i]] <- demonER(fit = qdpm2_quant[[i]], metric = "qDPM_p2", alpha = 0.1, quantile = quantiles[i])
  
}

r_mpd <- do.call(rbind, r_mpd)
r_qdpm0 <- do.call(rbind, r_qdpm0)
r_qdpm1 <- do.call(rbind, r_qdpm1)
r_qdpm2 <- do.call(rbind, r_qdpm2)

r_mpd_01 <- do.call(rbind, r_mpd_01)
r_qdpm0_01 <- do.call(rbind, r_qdpm0_01)
r_qdpm1_01 <- do.call(rbind, r_qdpm1_01)
r_qdpm2_01 <- do.call(rbind, r_qdpm2_01)

er_phylo_quantile <- rbind(r_mpd, r_qdpm0, r_qdpm1, r_qdpm2, 
                           r_mpd_01, r_qdpm0_01, r_qdpm1_01, r_qdpm2_01)

write.csv(er_phylo_quantile, file = "output/EviRatio/ER_phylo_continental_quantile.csv")

##### Load fitted models for the taxonomic dimension and run demonER #####
source("R/NEON_diversity/R/Functions/demonER.R")

quantiles <- c(0.05, 0.1, 0.25, 0.50, 0.75, 0.90, 0.95)
nQuant <- length(quantiles)
metrics <- c("Smsd", "Sss", "H", "D")

# Load data 
# SR ~ MSD
load("Results/Regressions/QuantileREG/Trait/quant_trait_SR_site_0595.RData")
srs_0595 <- quantREG_sr
load("Results/Regressions/QuantileREG/Trait/quant_trait_SR_site.RData")
srs_2575 <- quantREG_sr

msd_quant <- list(srs_0595[[1]], srs_0595[[2]], srs_2575[[1]], srs_2575[[2]], 
                  srs_2575[[3]], srs_0595[[3]], srs_0595[[4]])

load("Results/Regressions/QuantileREG/Taxonomy/quant_taxo_SR_site_0595.RData")
sr_0595 <- quantREG_sr
load("Results/Regressions/QuantileREG/Taxonomy/quant_taxo_SR_site.RData")
sr_2575 <- quantREG_sr

sr_quant <- list(sr_0595[[1]], sr_0595[[2]], sr_2575[[1]], sr_2575[[2]], 
                 sr_2575[[3]], sr_0595[[3]], sr_0595[[4]])

load("Results/Regressions/QuantileREG/Taxonomy/quant_taxo_H_site_0595.RData")
h_0595 <- quantREG_h
load("Results/Regressions/QuantileREG/Taxonomy/quant_taxo_H_site.RData")
h_2575 <- quantREG_h

h_quant <- list(h_0595[[1]], h_0595[[2]], h_2575[[1]], h_2575[[2]], 
                h_2575[[3]], h_0595[[3]], h_0595[[4]])

load("Results/Regressions/QuantileREG/Taxonomy/quant_taxo_D_site_0595.RData")
d_0595 <- quantREG_d
load("Results/Regressions/QuantileREG/Taxonomy/quant_taxo_D_site.RData")
d_2575 <- quantREG_d

d_quant <- list(d_0595[[1]], d_0595[[2]], d_2575[[1]], d_2575[[2]], 
                d_2575[[3]], d_0595[[3]], d_0595[[4]])

### Run ER for quantiles 
r_srs <- list()
r_sr <- list()
r_h <- list()
r_d <- list()

r_srs_01 <- list()
r_sr_01 <- list()
r_h_01 <- list()
r_d_01 <- list()

for(i in 1:length(quantiles)) {
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.05"))
  
  r_srs[[i]] <- demonER(fit = msd_quant[[i]], metric = "Smsd", alpha = 0.05, quantile = quantiles[i])
  r_sr[[i]] <- demonER(fit = sr_quant[[i]], metric = "S_q0", alpha = 0.05, quantile = quantiles[i])
  r_h[[i]] <- demonER(fit = h_quant[[i]], metric = "H_q1", alpha = 0.05, quantile = quantiles[i])
  r_d[[i]] <- demonER(fit = d_quant[[i]], metric = "D_q2", alpha = 0.05, quantile = quantiles[i])
  
  print(paste0("Run ER for quantile ", quantiles[i], " alpha = 0.1"))
  
  r_srs_01[[i]] <- demonER(fit = msd_quant[[i]], metric = "Smsd", alpha = 0.1, quantile = quantiles[i])
  r_sr_01[[i]] <- demonER(fit = sr_quant[[i]], metric = "S_q0", alpha = 0.1, quantile = quantiles[i])
  r_h_01[[i]] <- demonER(fit = h_quant[[i]], metric = "H_q1", alpha = 0.1, quantile = quantiles[i])
  r_d_01[[i]] <- demonER(fit = d_quant[[i]], metric = "D_q2", alpha = 0.1, quantile = quantiles[i])
  
}

r_srs <- do.call(rbind, r_srs)
r_sr <- do.call(rbind, r_sr)
r_h <- do.call(rbind, r_h)
r_d <- do.call(rbind, r_d)

r_srs_01 <- do.call(rbind, r_srs_01)
r_sr_01 <- do.call(rbind, r_sr_01)
r_h_01 <- do.call(rbind, r_h_01)
r_d_01 <- do.call(rbind, r_d_01)

er_taxo_quantile <- rbind(r_srs, r_sr, r_h, r_d, 
                           r_srs_01, r_sr_01, r_h_01, r_d_01)

write.csv(er_taxo_quantile, file = "output/EviRatio/ER_taxo_continental_quantile.csv")
