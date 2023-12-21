library(brms)
library(tidyverse)

########## ---------- Trait dimension ---------- ##########

##### Load data #####
## Q0
trait_spec_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q0.csv") 
## Q1
trait_spec_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q1.csv") 
## Q2
trait_spec_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q2.csv") 

##### Trait metric correlations ##### 

## MFD ~ MFDz
mfd <- brm(data = trait_spec_NEON_table_q0, 
           family = student, 
           bf(mvbind(MFD_trait, MFDz_trait) ~ 1) + set_rescor(TRUE), 
           iter = 2000, warmup = 500, chains = 4, cores = 4, 
           backend = "cmdstanr")

## MFDz ~ qDTM
mfdz_qdtm0 <- brm(data = trait_spec_NEON_table_q0, 
                 family = student, 
                 bf(mvbind(MFDz_trait, qDTM_trait) ~ 1) + set_rescor(TRUE), 
                 iter = 2000, warmup = 500, chains = 4, cores = 4, 
                 backend = "cmdstanr")

mfdz_qdtm1 <- brm(data = trait_spec_NEON_table_q1, 
                 family = student, 
                 bf(mvbind(MFDz_trait, qDTM_trait) ~ 1) + set_rescor(TRUE), 
                 iter = 2000, warmup = 500, chains = 4, cores = 4, 
                 backend = "cmdstanr")

mfdz_qdtm2 <- brm(data = trait_spec_NEON_table_q2, 
                 family = student, 
                 bf(mvbind(MFDz_trait, qDTM_trait) ~ 1) + set_rescor(TRUE), 
                 iter = 2000, warmup = 500, chains = 4, cores = 4, 
                 backend = "cmdstanr")

## MFD ~ qDTM
mfd_qdtm0 <- brm(data = trait_spec_NEON_table_q0, 
                  family = student, 
                  bf(mvbind(MFD_trait, qDTM_trait) ~ 1) + set_rescor(TRUE), 
                  iter = 2000, warmup = 500, chains = 4, cores = 4, 
                  backend = "cmdstanr")

mfd_qdtm1 <- brm(data = trait_spec_NEON_table_q1, 
                  family = student, 
                  bf(mvbind(MFD_trait, qDTM_trait) ~ 1) + set_rescor(TRUE), 
                  iter = 2000, warmup = 500, chains = 4, cores = 4, 
                  backend = "cmdstanr")

mfd_qdtm2 <- brm(data = trait_spec_NEON_table_q2, 
                  family = student, 
                  bf(mvbind(MFD_trait, qDTM_trait) ~ 1) + set_rescor(TRUE), 
                  iter = 2000, warmup = 500, chains = 4, cores = 4, 
                  backend = "cmdstanr")

traitcors <- list(mfd, mfd_qdtm0, mfd_qdtm1, mfd_qdtm2, 
                  mfdz_qdtm0, mfdz_qdtm1, mfdz_qdtm2)

save(traitcors, 
     file = "output/metricCor/Reanalyses/trait_cors.RData")

########## ---------- Phylogenetic dimension ---------- ########## 

##### Load data #####
## Q0
phylo_spec_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q0.csv")
## Q1
phylo_spec_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q1.csv")
## Q2
phylo_spec_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q2.csv")

##### Phylogenetic metric correlations #####

## MPD ~ MPDz
mpd <- brm(data = phylo_spec_NEON_table_q0, 
           family = student, 
           bf(mvbind(MPD_phylo, MPDz_phylo) ~ 1) + set_rescor(TRUE), 
           iter = 2000, warmup = 500, chains = 4, cores = 4, 
           backend = "cmdstanr")

## MPDz ~ qDTM
mpdz_qdtm0 <- brm(data = phylo_spec_NEON_table_q0, 
                  family = student, 
                  bf(mvbind(MPDz_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                  iter = 2000, warmup = 500, chains = 4, cores = 4, 
                  backend = "cmdstanr")

mpdz_qdtm1 <- brm(data = phylo_spec_NEON_table_q1, 
                  family = student, 
                  bf(mvbind(MPDz_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                  iter = 2000, warmup = 500, chains = 4, cores = 4, 
                  backend = "cmdstanr")

mpdz_qdtm2 <- brm(data = phylo_spec_NEON_table_q2, 
                  family = student, 
                  bf(mvbind(MPDz_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                  iter = 2000, warmup = 500, chains = 4, cores = 4, 
                  backend = "cmdstanr")

## MPD ~ qDTM
mpd_qdtm0 <- brm(data = phylo_spec_NEON_table_q0, 
                 family = student, 
                 bf(mvbind(MPD_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                 iter = 2000, warmup = 500, chains = 4, cores = 4, 
                 backend = "cmdstanr")

mpd_qdtm1 <- brm(data = phylo_spec_NEON_table_q1, 
                 family = student, 
                 bf(mvbind(MPD_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                 iter = 2000, warmup = 500, chains = 4, cores = 4, 
                 backend = "cmdstanr")

mpd_qdtm2 <- brm(data = phylo_spec_NEON_table_q2, 
                 family = student, 
                 bf(mvbind(MPD_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                 iter = 2000, warmup = 500, chains = 4, cores = 4, 
                 backend = "cmdstanr")

phylocors <- list(mpd, mpd_qdtm0, mpd_qdtm1, mpd_qdtm2, 
                  mpdz_qdtm0, mpdz_qdtm1, mpdz_qdtm2)

save(phylocors, 
     file = "output/metricCor/Reanalyses/phylo_cors.RData")

## PD ~ PDz
pd <- brm(data = phylo_spec_NEON_table_q0, 
           family = student, 
           bf(mvbind(PD_phylo, PDz_phylo) ~ 1) + set_rescor(TRUE), 
           iter = 2000, warmup = 500, chains = 4, cores = 4, 
           backend = "cmdstanr")

## PDz ~ qDTM
pdz_qdtm0 <- brm(data = phylo_spec_NEON_table_q0, 
                  family = student, 
                  bf(mvbind(PDz_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                  iter = 2000, warmup = 500, chains = 4, cores = 4, 
                  backend = "cmdstanr")

pdz_qdtm1 <- brm(data = phylo_spec_NEON_table_q1, 
                  family = student, 
                  bf(mvbind(PDz_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                  iter = 2000, warmup = 500, chains = 4, cores = 4, 
                  backend = "cmdstanr")

pdz_qdtm2 <- brm(data = phylo_spec_NEON_table_q2, 
                  family = student, 
                  bf(mvbind(PDz_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                  iter = 2000, warmup = 500, chains = 4, cores = 4, 
                  backend = "cmdstanr")

## PD ~ qDTM
pd_qdtm0 <- brm(data = phylo_spec_NEON_table_q0, 
                 family = student, 
                 bf(mvbind(PD_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                 iter = 2000, warmup = 500, chains = 4, cores = 4, 
                 backend = "cmdstanr")

pd_qdtm1 <- brm(data = phylo_spec_NEON_table_q1, 
                 family = student, 
                 bf(mvbind(PD_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                 iter = 2000, warmup = 500, chains = 4, cores = 4, 
                 backend = "cmdstanr")

pd_qdtm2 <- brm(data = phylo_spec_NEON_table_q2, 
                 family = student, 
                 bf(mvbind(PD_phylo, qDTM_phylo) ~ 1) + set_rescor(TRUE), 
                 iter = 2000, warmup = 500, chains = 4, cores = 4, 
                 backend = "cmdstanr")

phylocors_PD <- list(pd, pd_qdtm0, pd_qdtm1, pd_qdtm2, 
                     pdz_qdtm0, pdz_qdtm1, pdz_qdtm2)

save(phylocors_PD, 
     file = "output/metricCor/Reanalyses/phylo_cors_PD.RData")

########## ---------- Taxonomic dimension ---------- ########## 

##### Load data #####
taxo_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/taxo_spec_SS_NEON.csv") %>% 
  mutate(siteID = Site)

##### Taxonomic metric correlations #####
## S ~ H
s_h <- brm(data = taxo_spec_SS_NEON, 
           family = student, 
           bf(mvbind(SRich_ground, Shannon_ground) ~ 1) + set_rescor(TRUE), 
           iter = 2000, warmup = 500, chains = 4, cores = 4, 
           backend = "cmdstanr")

## S ~ D
s_d <- brm(data = taxo_spec_SS_NEON, 
           family = student, 
           bf(mvbind(SRich_ground, Simpson_ground) ~ 1) + set_rescor(TRUE), 
           iter = 2000, warmup = 500, chains = 4, cores = 4, 
           backend = "cmdstanr")

## H ~ D
h_d  <- brm(data = taxo_spec_SS_NEON, 
            family = student, 
            bf(mvbind(Shannon_ground, Simpson_ground) ~ 1) + set_rescor(TRUE), 
            iter = 2000, warmup = 500, chains = 4, cores = 4, 
            backend = "cmdstanr")

taxocors <- list(s_h, s_d, h_d)

save(taxocors, 
     file = "output/metricCor/Reanalyses/taxo_cors.RData")

########## ---------- Spectral dimension ---------- ########## 

##### Load data ##### 

### Distance based metrics
phylo_spec_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q0.csv") 

### Spectral species based metrics
taxo_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/taxo_spec_SS_NEON.csv") 

##### Spectral distance metric correlations #####

## MSD ~ qDSM
msd_qdsm <- brm(data = phylo_spec_NEON_table_q0, 
           family = student, 
           bf(mvbind(MSD_spec, qDTM_spec) ~ 1) + set_rescor(TRUE), 
           iter = 2000, warmup = 500, chains = 4, cores = 4, 
           backend = "cmdstanr")

## Additional correlations between metrics of spectral dimension 
## SD ~ qDSM
sd_qdsm <- brm(data = phylo_spec_NEON_table_q0, 
                family = student, 
                bf(mvbind(SD_spec, qDTM_spec) ~ 1) + set_rescor(TRUE), 
                iter = 2000, warmup = 500, chains = 4, cores = 4, 
                backend = "cmdstanr")

## SD ~ MSD
sd_msd <- brm(data = phylo_spec_NEON_table_q0, 
               family = student, 
               bf(mvbind(SD_spec, MSD_spec) ~ 1) + set_rescor(TRUE), 
               iter = 2000, warmup = 500, chains = 4, cores = 4, 
               backend = "cmdstanr")

speccors_dis <- list(sd_msd, sd_qdsm, msd_qdsm)

save(speccors_dis, 
     file = "output/metricCor/Reanalyses/spec_cors_DIS.RData")

##### Spectral species metric correlations #####

## S ~ H
s_h_spec <- brm(data = taxo_spec_SS_NEON, 
                family = student, 
                bf(mvbind(SRich_spec, Shannon_spec) ~ 1) + set_rescor(TRUE), 
                iter = 2000, warmup = 500, chains = 4, cores = 4, 
                backend = "cmdstanr")

## S ~ D
s_d_spec <- brm(data = taxo_spec_SS_NEON, 
               family = student, 
               bf(mvbind(SRich_spec, Simpson_spec) ~ 1) + set_rescor(TRUE), 
               iter = 2000, warmup = 500, chains = 4, cores = 4, 
               backend = "cmdstanr")

## H ~ D
h_d_spec <- brm(data = taxo_spec_SS_NEON, 
                family = student, 
                bf(mvbind(Shannon_spec, Simpson_spec) ~ 1) + set_rescor(TRUE), 
                iter = 2000, warmup = 500, chains = 4, cores = 4, 
                backend = "cmdstanr")

speccors_ss <- list(s_h_spec, s_d_spec, h_d_spec)

save(speccors_ss, 
     file = "output/metricCor/Reanalyses/spec_cors_SS.RData")

##### Get metric correlations #####
library(tidyverse)

### Auxiliary function
getCors <- function(fits, dimension, q = q) { 
  
  ### main library
  library(brms) 
  
  ### store results
  resCors <- list()
  
  for(i in 1:length(fits)) {
    
    ### Get correlations 
    res <- summary(fits[[i]])$rescor_pars
    res$metrics <- rownames(res)
    rownames(res) <- NULL
    
    resCors[[i]] <- res
    
  }
  
  ### combine results
  resCors <- do.call(rbind, resCors) 
  resCors$dimension <- dimension
  resCors$q <- q 
  
  resCors <- resCors %>% 
    select(dimension, metrics, q, everything())
  
  ### output
  return(resCors)
  
}

##### Load correlations #####

## phylogeny MPD
load("output/metricCor/Reanalyses/phylo_cors.RData") 
## phylogeny PD
load("output/metricCor/Reanalyses/phylo_cors_PD.RData") 
## trait MTD
load("output/metricCor/Reanalyses/trait_cors.RData") 
## Spectral distance
load("output/metricCor/Reanalyses/spec_cors_DIS.RData") 
## Spectral SS
load("output/metricCor/Reanalyses/spec_cors_SS.RData") 
## taxonomy
load("output/metricCor/Reanalyses/taxo_cors.RData") 

## Phylogeny MPD
corPhylo_MPD <- getCors(fits = phylocors, 
                        dimension = "phylogeny", 
                        q = c("q0", "q0", "q1", "q2", "q0", "q1", "q2")) 

## Phylogeny PD
corPhylo_PD <- getCors(fits = phylocors_PD, 
                        dimension = "phylogeny", 
                        q = c("q0", "q0", "q1", "q2", "q0", "q1", "q2")) 

## Trait MTD
corTrait <- getCors(fits = traitcors, 
                    dimension = "trait", 
                    q = c("q0", "q0", "q1", "q2", "q0", "q1", "q2"))

## Taxonomy 
corTaxo <- getCors(fits = taxocors, 
                   dimension = "taxonomy", 
                   q = c("q01", "q02", "q12"))

## Spectral distance
corSpec_dis <- getCors(fits = speccors_dis, 
                       dimension = "spectral", 
                       q = c("q0", "q0", "q0"))

## Spectral species 
corSpec_ss <- getCors(fits = speccors_ss, 
                       dimension = "spectral", 
                       q = c("q01", "q02", "q12"))

##### Combine results #####
corMetric <- bind_rows(corTaxo, 
                       corTrait, 
                       corPhylo_MPD, 
                       corPhylo_PD, 
                       corSpec_dis, 
                       corSpec_ss
)

### Save results
write_csv(corMetric, 
          file = "output/metricCor/Reanalyses/corMetrics.csv")
