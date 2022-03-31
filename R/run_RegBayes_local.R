##### Phylogenetic dimension #####
library(tidyverse)
source("R/NEON_diversity/R/Functions/GOD_Bayes_local.R")

load("Results/RegDATA/SAM/phylo_spec_DIS_SAM_NEON.RData")

### Inits
QS <- c("q0", "q1", "q2", "q3")
nQ <- length(QS) 
sites <- unique(phylo_spec_NEON_table_SAM_q0$Site) 
nSites <- length(sites)

res_phylo_q0 <- list()
res_phylo_q1 <- list()
res_phylo_q2 <- list()

for(i in 1:nSites) { 
  print(sites[i])
  
  ### Q0
  dataPhyQ0 <- phylo_spec_NEON_table_SAM_q0 %>% 
    filter(Site == sites[i]) %>% 
    select(Site, plotID, PD, PDz, MPD, MPDz, M_p, mPrime_p, qDTM_p, 
           SD, SD2, MSD, MSD2, M_s, mPrime_s, qDTM_s)
  
  reg_q0 <- god_BayReg_local_phylo(resMetrics = dataPhyQ0, 
                                       Q = QS[1], nMetrics = 7, site = sites[i], 
                                       scale = TRUE, dimension = "phylogeny", 
                                       nChains = 4, nIters = 2000, nCores = 24, 
                                       engine = "cmdstanr", 
                                       control = list(adapt_delta = 0.99) 
                                   ) 
  res_phylo_q0[[i]] <- reg_q0
  
  ### Q1
  dataPhyQ1 <- phylo_spec_NEON_table_SAM_q1 %>% 
    filter(Site == sites[i]) %>% 
    select(Site, plotID, PD, PDz, MPD, MPDz, M_p, mPrime_p, qDTM_p, 
           SD, SD2, MSD, MSD2, M_s, mPrime_s, qDTM_s)
  
  reg_q1 <- god_BayReg_local_phylo(resMetrics = dataPhyQ1, 
                                       Q = QS[2], nMetrics = 7, site = sites[i], 
                                       scale = TRUE, dimension = "phylogeny", 
                                       nChains = 4, nIters = 2000, nCores = 24, 
                                       engine = "cmdstanr", 
                                       control = list(adapt_delta = 0.99) 
                                   ) 
  res_phylo_q1[[i]] <- reg_q1
  
  ### Q2
  dataPhyQ2 <- phylo_spec_NEON_table_SAM_q2 %>% 
    filter(Site == sites[i]) %>% 
    select(Site, plotID, PD, PDz, MPD, MPDz, M_p, mPrime_p, qDTM_p, 
           SD, SD2, MSD, MSD2, M_s, mPrime_s, qDTM_s)
  
  reg_q2 <- god_BayReg_local_phylo(resMetrics = dataPhyQ2, 
                                       Q = QS[3], nMetrics = 7, site = sites[i], 
                                       scale = TRUE, dimension = "phylogeny", 
                                       nChains = 4, nIters = 2000, nCores = 24, 
                                       engine = "cmdstanr", 
                                       control = list(adapt_delta = 0.99) 
                                   ) 
  res_phylo_q2[[i]] <- reg_q2
  
  gc()
}

dir.create("Results/Regressions/phylo-spec/site")

save(res_phylo_q0, res_phylo_q1, res_phylo_q2, 
     file = "Results/Regressions/phylo-spec/site/output_phylo_site.RData")

rm(list = ls())

##### Trait dimension #####
library(tidyverse)
source("R/NEON_diversity/R/Functions/GOD_Bayes_local.R")

load("Results/RegDATA/SAM/trait_spec_DIS_SAM_NEON.RData")

### Inits
QS <- c("q0", "q1", "q2", "q3")
nQ <- length(QS) 
sites <- unique(trait_spec_NEON_table_SAM_q0$Site) 
nSites <- length(sites)

res_trait_dis_q0 <- list()
res_trait_dis_q1 <- list()
res_trait_dis_q2 <- list()

for(i in 1:nSites) { 
  
  print(sites[i])
  
  ### Q0
  data_q0 <- trait_spec_NEON_table_SAM_q0 %>% 
    filter(Site == sites[i]) %>% 
    select(Site, plotID, MFD, MFDz, SR, M_f, mPrime_f, qDTM_f, 
           MSD, MSD2, MSm, M_s, mPrime_s, qDTM_s)
  
  reg_dis_q0 <- god_BayReg_trait_local_dis(resMetrics = data_q0, 
                                           Q = QS[1], nMetrics = 6, site = sites[i], 
                                           scale = TRUE, dimension = "trait", 
                                           nChains = 4, nIters = 2000, 
                                           nCores = 24, engine = "cmdstanr") 
  res_trait_dis_q0[[i]] <- reg_dis_q0 
  
  ### Q1
  data_q1 <- trait_spec_NEON_table_SAM_q1 %>% 
    filter(Site == sites[i]) %>% 
    select(Site, plotID, MFD, MFDz, SR, M_f, mPrime_f, qDTM_f, 
           MSD, MSD2, MSm, M_s, mPrime_s, qDTM_s)
  
  reg_dis_q1 <- god_BayReg_trait_local_dis(resMetrics = data_q1, 
                                           Q = QS[2], nMetrics = 6, site = sites[i], 
                                           scale = TRUE, dimension = "trait", 
                                           nChains = 4, nIters = 2000, 
                                           nCores = 24, engine = "cmdstanr") 
  res_trait_dis_q1[[i]] <- reg_dis_q1 
  
  ### Q2
  data_q2 <- trait_spec_NEON_table_SAM_q2 %>% 
    filter(Site == sites[i]) %>% 
    select(Site, plotID, MFD, MFDz, SR, M_f, mPrime_f, qDTM_f, 
           MSD, MSD2, MSm, M_s, mPrime_s, qDTM_s)
  
  reg_dis_q2 <- god_BayReg_trait_local_dis(resMetrics = data_q2, 
                                           Q = QS[3], nMetrics = 6, site = sites[i], 
                                           scale = TRUE, dimension = "trait", 
                                           nChains = 4, nIters = 2000, 
                                           nCores = 24, engine = "cmdstanr") 
  res_trait_dis_q2[[i]] <- reg_dis_q2 
  
  gc()
}

dir.create("Results/Regressions/trait-spec/site")

save(res_trait_dis_q0, res_trait_dis_q1, res_trait_dis_q2, 
     file = "Results/Regressions/trait-spec/site/output_trait_dis_site.RData")

rm(list = ls())

##### Trait dimension #####
source("R/NEON_diversity/R/Functions/GOD_Bayes_local.R")

load("Results/RegDATA/SAM/trait_spec_MET_SAM_NEON.RData")

### Inits
sites <- unique(trait_spec_NEON_table_MET_SAM_q1$Site) 

nSites <- length(sites)

res_trait_met <- list()

for(i in 1:nSites) {
  
  print(sites[i]) 
  
  data_met <- trait_spec_NEON_table_MET_SAM_q1 %>% 
    filter(Site == sites[i])
  
  reg_met <- god_BayReg_trait_local_met(resMetrics = data_met, 
                                       Q = "q13", nMetrics = 5, site = sites[i], 
                                       scale = TRUE, dimension = "trait", 
                                       nChains = 4, nIters = 2000, 
                                       nCores = 24, engine = "cmdstanr") 
  res_trait_met[[i]] <- reg_met
  
}

save(res_trait_met, 
     file = "Results/Regressions/trait-spec/site/output_trait_met_site.RData")

rm(list = ls())

gc()

##### Taxonomic dimension #####
source("R/NEON_diversity/R/Functions/GOD_Bayes_local.R")

load("Results/RegDATA/taxo_spec_ALPHA_NEON.RData")

### Inits
sites <- unique(taxo_spec_thresh_NEON_table$Site) 

nSites <- length(sites)

res_taxo <- list()

for(i in 1:nSites) { 
  
  print(sites[i]) 
  
  data_taxo <- taxo_spec_thresh_NEON_table %>% 
    filter(Site == sites[i]) 
  
  reg_tx <- god_BayReg_alpha_local_taxo(resMetrics = data_taxo, 
                                        nMetrics = 13, scale = FALSE, 
                                        site = sites[i], dimension = "taxonomy", 
                                        nChains = 4, nIters = 2000, 
                                        nCores = 24, engine = "cmdstanr") 
  res_taxo[[i]] <- reg_tx
  
}

save(res_taxo, 
     file = "Results/Regressions/taxo-spec/site/output_taxo_site.RData")
