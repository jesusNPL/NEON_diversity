
##### ------------------ Phylogenetic dimension R2 -------------------- #####
##### Extract data for plotting #####

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirPhy <- "Results/Regressions/phylo-spec/SAM/Scale/"

fits_phy <- list.files(dirPhy)

Qs <- c("q0", "q1", "q2", "q3")

r2_lst <- list()

for(i in 1:length(fits_phy)) {
  
  fit_phy <- fits_phy[i]
  
  load(paste0(dirPhy, fit_phy))
  
  r2 <- getR2Estimates(fits = res, robust = TRUE, 
                       dimension = "phylogeny", Q = Qs[i], 
                       level = "all") 
  r2_lst[[i]] <- r2
  
}

r2_table_phy <- do.call(rbind, r2_lst)

##### By habitat #####
list.files("Results/Regressions/phylo-spec/Habitat/Scale/")

dirPhy_habitat <- "Results/Regressions/phylo-spec/Habitat/Scale/"

prefix <- "reg_Phylo_Spec_DIS_"

Qs <- c("q0", "q1", "q2", "q3")

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

r2_q0_lst <- list()
r2_q1_lst <- list()
r2_q2_lst <- list()
r2_q3_lst <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirPhy_habitat, prefix, Qs[1], "_", habitat[i], ".RData")) 
  
  r2q0 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "phylogeny", Q = Qs[1], 
                         level = habitat[i]) 
  r2_q0_lst[[i]] <- r2q0 
  
  # q1
  load(paste0(dirPhy_habitat, prefix, Qs[2], "_", habitat[i], ".RData")) 
  
  r2q1 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "phylogeny", Q = Qs[2], 
                         level = habitat[i]) 
  r2_q1_lst[[i]] <- r2q1
  
  # q2
  load(paste0(dirPhy_habitat, prefix, Qs[3], "_", habitat[i], ".RData")) 
  
  r2q2 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "phylogeny", Q = Qs[3], 
                         level = habitat[i]) 
  r2_q2_lst[[i]] <- r2q2
  
  # q3
  load(paste0(dirPhy_habitat, prefix, Qs[4], "_", habitat[i], ".RData")) 
  
  r2q3 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "phylogeny", Q = Qs[4], 
                         level = habitat[i]) 
  r2_q3_lst[[i]] <- r2q3
  
}

r2_q0_table <- do.call(rbind, r2_q0_lst)
r2_q1_table <- do.call(rbind, r2_q1_lst)
r2_q2_table <- do.call(rbind, r2_q2_lst)
r2_q3_table <- do.call(rbind, r2_q3_lst)

##### Combine all tables #####
library(tidyverse)

r2_table_phylo <- rbind(r2_table_phy, r2_q0_table, r2_q1_table, 
                        r2_q2_table, r2_q3_table)

write.csv(r2_table_phylo, 
          file = "Results/Regressions/phylo-spec/Estimations_R2_phylo.csv")


##### ------------------ Trait dimension R2 -------------------- #####
##### Extract data for plotting #####

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirTrt <- "Results/Regressions/trait-spec/SAM/Scale/"

fits_trt <- list.files(dirTrt)

Qs <- c("q0", "q1", "q2", "q3")

r2_lst_dis <- list()
r2_lst_met <- list()

for(i in 1:length(Qs)) {
  # Distance based metrics
  print(fits_trt[i])
  
  fit_trt_dis <- fits_trt[i]
  
  load(paste0(dirTrt, fit_trt_dis))
  
  r2dis <- getR2Estimates(fits = res, robust = TRUE, 
                       dimension = "trait", Q = Qs[i], 
                       level = "all") 
  r2_lst_dis[[i]] <- r2dis
  
  # Classic metrics
  print(fits_trt[i + 4])
  
  fit_trt_met <- fits_trt[i + 4]
  
  load(paste0(dirTrt, fit_trt_met))
  
  r2met <- getR2Estimates(fits = res, robust = TRUE, 
                          dimension = "trait", Q = Qs[i], 
                          level = "all") 
  r2_lst_met[[i]] <- r2met
  
}

r2_table_trt_dis <- do.call(rbind, r2_lst_dis)
r2_table_trt_met <- do.call(rbind, r2_lst_met)

##### By habitat #####
list.files("Results/Regressions/trait-spec/Habitat/DIS/")

dirTrt_habitat_dis <- "Results/Regressions/trait-spec/Habitat/DIS/"

prefixDIS <- "reg_Trait_Spec_DIS_SAM_scaled_"

Qs <- c("q0", "q1", "q2", "q3")

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

### Distance based metrics
r2_q0_lst_dis <- list()
r2_q1_lst_dis <- list()
r2_q2_lst_dis <- list()
r2_q3_lst_dis <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[1], "_", habitat[i], ".RData")) 
  
  r2q0 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "trait", Q = Qs[1], 
                         level = habitat[i]) 
  r2_q0_lst_dis[[i]] <- r2q0 
  
  # q1
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[2], "_", habitat[i], ".RData")) 
  
  r2q1 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "trait", Q = Qs[2], 
                         level = habitat[i]) 
  r2_q1_lst_dis[[i]] <- r2q1
  
  # q2
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[3], "_", habitat[i], ".RData")) 
  
  r2q2 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "trait", Q = Qs[3], 
                         level = habitat[i]) 
  r2_q2_lst_dis[[i]] <- r2q2
  
  # q3
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[4], "_", habitat[i], ".RData")) 
  
  r2q3 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "trait", Q = Qs[4], 
                         level = habitat[i]) 
  r2_q3_lst_dis[[i]] <- r2q3
  
}

r2_q0_table_dis <- do.call(rbind, r2_q0_lst_dis)
r2_q1_table_dis <- do.call(rbind, r2_q1_lst_dis)
r2_q2_table_dis <- do.call(rbind, r2_q2_lst_dis)
r2_q3_table_dis <- do.call(rbind, r2_q3_lst_dis)

#### Classic metrics  
list.files("Results/Regressions/trait-spec/Habitat/MET/")

dirTrt_habitat_met <- "Results/Regressions/trait-spec/Habitat/MET/"

prefixMET <- "reg_Trait_Spec_MET_SAM_scaled_"

r2_q0_lst_met <- list()
r2_q1_lst_met <- list()
r2_q2_lst_met <- list()
r2_q3_lst_met <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[1], "_", habitat[i], ".RData")) 
  
  r2q0 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "trait", Q = Qs[1], 
                         level = habitat[i]) 
  r2_q0_lst_met[[i]] <- r2q0 
  
  # q1
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[2], "_", habitat[i], ".RData")) 
  
  r2q1 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "trait", Q = Qs[2], 
                         level = habitat[i]) 
  r2_q1_lst_met[[i]] <- r2q1
  
  # q2
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[3], "_", habitat[i], ".RData")) 
  
  r2q2 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "trait", Q = Qs[3], 
                         level = habitat[i]) 
  r2_q2_lst_met[[i]] <- r2q2
  
  # q3
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[4], "_", habitat[i], ".RData")) 
  
  r2q3 <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "trait", Q = Qs[4], 
                         level = habitat[i]) 
  r2_q3_lst_met[[i]] <- r2q3
  
}

r2_q0_table_met <- do.call(rbind, r2_q0_lst_met)
r2_q1_table_met <- do.call(rbind, r2_q1_lst_met)
r2_q2_table_met <- do.call(rbind, r2_q2_lst_met)
r2_q3_table_met <- do.call(rbind, r2_q3_lst_met)

##### Combine all tables #####
library(tidyverse)

r2_table_trait_dis <- rbind(r2_table_trt_dis, r2_q0_table_dis, 
                            r2_q1_table_dis, r2_q2_table_dis, 
                            r2_q3_table_dis)

r2_table_trait_met <- rbind(r2_table_trt_met, r2_q0_table_met, 
                            r2_q1_table_met, r2_q2_table_met, 
                            r2_q3_table_met)

write.csv(r2_table_trait_dis, 
          file = "Results/Regressions/trait-spec/Estimations_R2_trait_dis.csv")

write.csv(r2_table_trait_met, 
          file = "Results/Regressions/trait-spec/Estimations_R2_trait_met.csv")

##### ------------------ Taxonomic dimension R2 -------------------- #####
##### Extract data for plotting #####

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirTax <- "Results/Regressions/taxo-spec/"

fits_tax <- list.files(dirTax)


load(paste0(dirPhy, fits_tax[22]))
  
r2_tax <- getR2Estimates(fits = res, robust = TRUE, 
                       dimension = "taxonomy", Q = "q13", 
                       level = "all") 
  
##### By habitat #####
list.files("Results/Regressions/taxo-spec/Habitat/")

dirTax_habitat <- "Results/Regressions/taxo-spec/Habitat/"

prefix <- "reg_Taxo_Spec_Alpha_THRESH_"

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

r2_tax_lst <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirTax_habitat, prefix, habitat[i], ".RData")) 
  
  r2tax <- getR2Estimates(fits = res, robust = TRUE, 
                         dimension = "taxonomy", Q = "q13", 
                         level = habitat[i]) 
  r2_tax_lst[[i]] <- r2tax 
  
}

r2_tax_table <- do.call(rbind, r2_tax_lst)

##### Combine all tables #####
r2_table_taxo <- rbind(r2_tax, r2_tax_table)

write.csv(r2_table_taxo, 
          file = "Results/Regressions/taxo-spec/Estimations_R2_taxo.csv")

##### --------------- Extract Betas PHYLOGENY ------------------ #####

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirPhy <- "Results/Regressions/phylo-spec/SAM/Scale/"

fits_phy <- list.files(dirPhy)

Qs <- c("q0", "q1", "q2", "q3")

beta_lst <- list()

for(i in 1:length(fits_phy)) {
  
  fit_phy <- fits_phy[i]
  
  load(paste0(dirPhy, fit_phy))
  
  betas <- getFixefPosterior(fits = res$fits, 
                             estimates = res$R2_robust, 
                             nDraws = 1000, 
                             dimension = "phylogeny", 
                             Q = Qs[i], 
                             level = "all") 
  beta_lst[[i]] <- betas
  
}

beta_table_phylo <- do.call(rbind, beta_lst)

##### By habitat #####
list.files("Results/Regressions/phylo-spec/Habitat/Scale/")

dirPhy_habitat <- "Results/Regressions/phylo-spec/Habitat/Scale/"

prefix <- "reg_Phylo_Spec_DIS_"

Qs <- c("q0", "q1", "q2", "q3")

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

beta_q0_lst <- list()
beta_q1_lst <- list()
beta_q2_lst <- list()
beta_q3_lst <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirPhy_habitat, prefix, Qs[1], "_", habitat[i], ".RData")) 
  
  betaq0 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
                              dimension = "phylogeny", 
                              Q = Qs[1], 
                              level = habitat[i]) 
  beta_q0_lst[[i]] <- betaq0 
  
  # q1
  load(paste0(dirPhy_habitat, prefix, Qs[2], "_", habitat[i], ".RData")) 
  
  betaq1 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
                              dimension = "phylogeny",
                              Q = Qs[2], 
                              level = habitat[i]) 
  beta_q1_lst[[i]] <- betaq1
  
  # q2
  load(paste0(dirPhy_habitat, prefix, Qs[3], "_", habitat[i], ".RData")) 
  
  betaq2 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
                              dimension = "phylogeny",
                              Q = Qs[3], 
                              level = habitat[i]) 
  beta_q2_lst[[i]] <- betaq2
  
  # q3
  load(paste0(dirPhy_habitat, prefix, Qs[4], "_", habitat[i], ".RData")) 
  
  betaq3 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
                              dimension = "phylogeny",
                              Q = Qs[4], 
                              level = habitat[i]) 
  beta_q3_lst[[i]] <- betaq3
  
}

beta_q0_table <- do.call(rbind, beta_q0_lst)
beta_q1_table <- do.call(rbind, beta_q1_lst)
beta_q2_table <- do.call(rbind, beta_q2_lst)
beta_q3_table <- do.call(rbind, beta_q3_lst)

##### Combine all tables #####
library(tidyverse)

beta_table_phylogeny <- rbind(beta_table_phylo, beta_q0_table, beta_q1_table, 
                              beta_q2_table, beta_q3_table)

write.csv(beta_table_phylogeny, 
          file = "Results/Regressions/phylo-spec/Estimations_BETA_phylo.csv")

##### --------------- Extract Betas TRAIT ------------------ #####

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirTrt <- "Results/Regressions/trait-spec/SAM/Scale/"

fits_trt <- list.files(dirTrt)

Qs <- c("q0", "q1", "q2", "q3")

beta_lst_dis <- list()
beta_lst_met <- list()

for(i in 1:length(Qs)) {
  
  # Distance based metrics
  print(fits_trt[i])

  fit_trt_dis <- fits_trt[i]
  
  load(paste0(dirTrt, fit_trt_dis))
  
  betas_dis <- getFixefPosterior(fits = res$fits, 
                             estimates = res$R2_robust, 
                             nDraws = 5000, 
                             dimension = "trait", 
                             Q = Qs[i], 
                             level = "all") 
  beta_lst_dis[[i]] <- betas_dis
  
  # classic metrics
  print(fits_trt[i + 4])
  
  fit_trt_met <- fits_trt[i + 4]
  
  load(paste0(dirTrt, fit_trt_met))
  
  betas_met <- getFixefPosterior(fits = res$fits, 
                                 estimates = res$R2_robust, 
                                 nDraws = 5000, 
                                 dimension = "trait", 
                                 Q = Qs[i], 
                                 level = "all") 
  beta_lst_met[[i]] <- betas_met
}

beta_table_trait_dis <- do.call(rbind, beta_lst_dis)
beta_table_trait_met <- do.call(rbind, beta_lst_met)

##### By habitat distance based #####

list.files("Results/Regressions/trait-spec/Habitat/DIS/")

dirTrt_habitat_dis <- "Results/Regressions/trait-spec/Habitat/DIS/"

prefix_dis <- "reg_Trait_Spec_DIS_SAM_scaled_"

Qs <- c("q0", "q1", "q2", "q3")

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

beta_q0_lst_dis <- list()
beta_q1_lst_dis <- list()
beta_q2_lst_dis <- list()
beta_q3_lst_dis <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirTrt_habitat_dis, prefix_dis, Qs[1], "_", habitat[i], ".RData")) 
  
  betaq0 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 5000, 
                              dimension = "trait", 
                              Q = Qs[1], 
                              level = habitat[i]) 
  beta_q0_lst_dis[[i]] <- betaq0 
  
  # q1
  load(paste0(dirTrt_habitat_dis, prefix_dis, Qs[2], "_", habitat[i], ".RData")) 
  
  betaq1 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 5000, 
                              dimension = "trait",
                              Q = Qs[2], 
                              level = habitat[i]) 
  beta_q1_lst_dis[[i]] <- betaq1
  
  # q2
  load(paste0(dirTrt_habitat_dis, prefix_dis, Qs[3], "_", habitat[i], ".RData")) 
  
  betaq2 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 5000, 
                              dimension = "trait",
                              Q = Qs[3], 
                              level = habitat[i]) 
  beta_q2_lst_dis[[i]] <- betaq2
  
  # q3
  load(paste0(dirTrt_habitat_dis, prefix_dis, Qs[4], "_", habitat[i], ".RData")) 
  
  betaq3 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 5000, 
                              dimension = "trait",
                              Q = Qs[4], 
                              level = habitat[i]) 
  beta_q3_lst_dis[[i]] <- betaq3
  
}

beta_q0_table_dis <- do.call(rbind, beta_q0_lst_dis)
beta_q1_table_dis <- do.call(rbind, beta_q1_lst_dis)
beta_q2_table_dis <- do.call(rbind, beta_q2_lst_dis)
beta_q3_table_dis <- do.call(rbind, beta_q3_lst_dis)

##### Combine all tables #####
beta_table_trait_dis_ALL <- rbind(beta_table_trait_dis, beta_q0_table_dis, 
                              beta_q1_table_dis, beta_q2_table_dis, beta_q3_table_dis)

write.csv(beta_table_trait_dis_ALL, 
          file = "Results/Regressions/trait-spec/Estimations_BETA_trait_dis.csv")

##### By habitat classic metrics #####

list.files("Results/Regressions/trait-spec/Habitat/MET/")

dirTrt_habitat_met <- "Results/Regressions/trait-spec/Habitat/MET/"

prefix_met <- "reg_Trait_Spec_MET_SAM_scaled_"

Qs <- c("q0", "q1", "q2", "q3")

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

beta_q0_lst_met <- list()
beta_q1_lst_met <- list()
beta_q2_lst_met <- list()
beta_q3_lst_met <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirTrt_habitat_met, prefix_met, Qs[1], "_", habitat[i], ".RData")) 
  
  betaq0 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 5000, 
                              dimension = "trait", 
                              Q = Qs[1], 
                              level = habitat[i]) 
  beta_q0_lst_met[[i]] <- betaq0 
  
  # q1
  load(paste0(dirTrt_habitat_met, prefix_met, Qs[2], "_", habitat[i], ".RData")) 
  
  betaq1 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 5000, 
                              dimension = "trait",
                              Q = Qs[2], 
                              level = habitat[i]) 
  beta_q1_lst_met[[i]] <- betaq1
  
  # q2
  load(paste0(dirTrt_habitat_met, prefix_met, Qs[3], "_", habitat[i], ".RData")) 
  
  betaq2 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 5000, 
                              dimension = "trait",
                              Q = Qs[3], 
                              level = habitat[i]) 
  beta_q2_lst_met[[i]] <- betaq2
  
  # q3
  load(paste0(dirTrt_habitat_met, prefix_met, Qs[4], "_", habitat[i], ".RData")) 
  
  betaq3 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 5000, 
                              dimension = "trait",
                              Q = Qs[4], 
                              level = habitat[i]) 
  beta_q3_lst_met[[i]] <- betaq3
  
}

beta_q0_table_met <- do.call(rbind, beta_q0_lst_met)
beta_q1_table_met <- do.call(rbind, beta_q1_lst_met)
beta_q2_table_met <- do.call(rbind, beta_q2_lst_met)
beta_q3_table_met <- do.call(rbind, beta_q3_lst_met)

##### Combine all tables #####
beta_table_trait_met_ALL <- rbind(beta_table_trait_met, beta_q0_table_met, 
                                  beta_q1_table_met, beta_q2_table_met, beta_q3_table_met)

write.csv(beta_table_trait_met_ALL, 
          file = "Results/Regressions/trait-spec/Estimations_BETA_trait_met.csv")

############ ------------ MAKE HYPOTHESIS --------------- ##############

##### Phylogenetic dimension #####
source("R/NEON_diversity/R/Functions/getBEffects.R")

dirPhy <- "Results/Regressions/phylo-spec/SAM/Scale/"

fits_phy <- list.files(dirPhy)

Qs <- c("q0", "q1", "q2", "q3")

h_lst <- list()

for(i in 1:length(fits_phy)) {
  
  fit_phy <- fits_phy[i]
  
  load(paste0(dirPhy, fit_phy))
  
  h <- makeHypothesis(fits = res$fits, 
                      estimates = res$R2_robust, 
                      dimension = "phylogeny", 
                      Q = Qs[i], 
                      level = "all") 
  h_lst[[i]] <- h
  
}

h_table_phy <- do.call(rbind, h_lst)

write.csv(h_table_phy, file = "output/hypotesis_phylo_allSites.csv")

##### By habitat #####
list.files("Results/Regressions/phylo-spec/Habitat/Scale/")

dirPhy_habitat <- "Results/Regressions/phylo-spec/Habitat/Scale/"

prefix <- "reg_Phylo_Spec_DIS_"

Qs <- c("q0", "q1", "q2", "q3")

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

h_q0_lst <- list()
h_q1_lst <- list()
h_q2_lst <- list()
h_q3_lst <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirPhy_habitat, prefix, Qs[1], "_", habitat[i], ".RData")) 
  
  hq0 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust, 
                        dimension = "phylogeny", Q = Qs[1], 
                        level = habitat[i]) 
  h_q0_lst[[i]] <- hq0 
  
  # q1
  load(paste0(dirPhy_habitat, prefix, Qs[2], "_", habitat[i], ".RData")) 
  
  hq1 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust,
                        dimension = "phylogeny", Q = Qs[2], 
                        level = habitat[i]) 
  h_q1_lst[[i]] <- hq1
  
  # q2
  load(paste0(dirPhy_habitat, prefix, Qs[3], "_", habitat[i], ".RData")) 
  
  hq2 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust,
                        dimension = "phylogeny", Q = Qs[3], 
                        level = habitat[i]) 
  h_q2_lst[[i]] <- hq2
  
  # q3
  load(paste0(dirPhy_habitat, prefix, Qs[4], "_", habitat[i], ".RData")) 
  
  hq3 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust,
                        dimension = "phylogeny", Q = Qs[4], 
                        level = habitat[i]) 
  h_q3_lst[[i]] <- hq3
  
}

h_q0_table <- do.call(rbind, h_q0_lst)
h_q1_table <- do.call(rbind, h_q1_lst)
h_q2_table <- do.call(rbind, h_q2_lst)
h_q3_table <- do.call(rbind, h_q3_lst)

##### Combine all tables #####
h_table_phylo_habitat <- rbind(h_q0_table, h_q1_table, 
                        h_q2_table, h_q3_table)

write.csv(h_table_phylo_habitat, 
          file = "output/hypotesis_phylo_habitat.csv")

##### Trait dimension #####

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirTrt <- "Results/Regressions/trait-spec/SAM/Scale/"

fits_trt <- list.files(dirTrt)

Qs <- c("q0", "q1", "q2", "q3")

h_lst_dis <- list()
h_lst_met <- list()

for(i in 1:length(Qs)) {
  # Distance based metrics
  print(fits_trt[i])
  
  fit_trt_dis <- fits_trt[i]
  
  load(paste0(dirTrt, fit_trt_dis))
  
  hdis <- makeHypothesis(fits = res$fits, 
                         estimates = res$R2_robust, 
                         dimension = "trait", Q = Qs[i], 
                         level = "all") 
  h_lst_dis[[i]] <- hdis
  
  # Classic metrics
  print(fits_trt[i + 4])
  
  fit_trt_met <- fits_trt[i + 4]
  
  load(paste0(dirTrt, fit_trt_met))
  
  hmet <- makeHypothesis(fits = res$fits, 
                         estimates = res$R2_robust, 
                         dimension = "trait", Q = Qs[i], 
                         level = "all") 
  h_lst_met[[i]] <- hmet
  
}

h_table_trt_dis <- do.call(rbind, h_lst_dis)
h_table_trt_met <- do.call(rbind, h_lst_met)

h_table_trait <- rbind(h_table_trt_dis, h_table_trt_met)

write.csv(h_table_trait, file = "output/hypotesis_trait_allSites.csv")

##### By habitat #####
list.files("Results/Regressions/trait-spec/Habitat/DIS/")

dirTrt_habitat_dis <- "Results/Regressions/trait-spec/Habitat/DIS/"

prefixDIS <- "reg_Trait_Spec_DIS_SAM_scaled_"

Qs <- c("q0", "q1", "q2", "q3")

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

### Distance based metrics
h_q0_lst_dis <- list()
h_q1_lst_dis <- list()
h_q2_lst_dis <- list()
h_q3_lst_dis <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[1], "_", habitat[i], ".RData")) 
  
  hq0 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust, 
                         dimension = "trait", Q = Qs[1], 
                         level = habitat[i]) 
  h_q0_lst_dis[[i]] <- hq0 
  
  # q1
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[2], "_", habitat[i], ".RData")) 
  
  hq1 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust, 
                        dimension = "trait", Q = Qs[2], 
                        level = habitat[i]) 
  h_q1_lst_dis[[i]] <- hq1
  
  # q2
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[3], "_", habitat[i], ".RData")) 
  
  hq2 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust, 
                        dimension = "trait", Q = Qs[3], 
                        level = habitat[i]) 
  h_q2_lst_dis[[i]] <- hq2
  
  # q3
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[4], "_", habitat[i], ".RData")) 
  
  hq3 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust, 
                        dimension = "trait", Q = Qs[4], 
                        level = habitat[i]) 
  h_q3_lst_dis[[i]] <- hq3
  
}

h_q0_table_dis <- do.call(rbind, h_q0_lst_dis)
h_q1_table_dis <- do.call(rbind, h_q1_lst_dis)
h_q2_table_dis <- do.call(rbind, h_q2_lst_dis)
h_q3_table_dis <- do.call(rbind, h_q3_lst_dis)

#### Classic metrics  
list.files("Results/Regressions/trait-spec/Habitat/MET/")

dirTrt_habitat_met <- "Results/Regressions/trait-spec/Habitat/MET/"

prefixMET <- "reg_Trait_Spec_MET_SAM_scaled_"

h_q0_lst_met <- list()
h_q1_lst_met <- list()
h_q2_lst_met <- list()
h_q3_lst_met <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[1], "_", habitat[i], ".RData")) 
  
  hq0 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust, 
                        dimension = "trait", Q = Qs[1], 
                        level = habitat[i]) 
  h_q0_lst_met[[i]] <- hq0 
  
  # q1
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[2], "_", habitat[i], ".RData")) 
  
  hq1 <- makeHypothesis(fits = res$fits, 
                         estimates = res$R2_robust, 
                         dimension = "trait", Q = Qs[2], 
                         level = habitat[i]) 
  h_q1_lst_met[[i]] <- hq1
  
  # q2
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[3], "_", habitat[i], ".RData")) 
  
  hq2 <- makeHypothesis(fits = res$fits, 
                        estimates = res$R2_robust, 
                        dimension = "trait", Q = Qs[3], 
                        level = habitat[i]) 
  h_q2_lst_met[[i]] <- hq2
  
  # q3
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[4], "_", habitat[i], ".RData")) 
  
  hq3 <- makeHypothesis(fits = res$fits, 
                         estimates = res$R2_robust, 
                         dimension = "trait", Q = Qs[4], 
                         level = habitat[i]) 
  h_q3_lst_met[[i]] <- hq3
  
}

h_q0_table_met <- do.call(rbind, h_q0_lst_met)
h_q1_table_met <- do.call(rbind, h_q1_lst_met)
h_q2_table_met <- do.call(rbind, h_q2_lst_met)
h_q3_table_met <- do.call(rbind, h_q3_lst_met)

##### Combine all tables #####

h_table_trait_dis <- rbind(h_q0_table_dis, h_q1_table_dis, 
                           h_q2_table_dis, h_q3_table_dis)

h_table_trait_met <- rbind(h_q0_table_met, h_q1_table_met, 
                           h_q2_table_met, h_q3_table_met)

write.csv(h_table_trait_dis, 
          file = "output/hypotesis_trait_habitat_dis.csv")

write.csv(h_table_trait_met, 
          file = "output/hypotesis_trait_habitat_met.csv")

############### ------------ GET FIXEF PHYLOGENY  ------------- ###############
source("R/NEON_diversity/R/Functions/getBEffects.R")

dirPhy <- "Results/Regressions/phylo-spec/SAM/Scale/"

fits_phy <- list.files(dirPhy)

Qs <- c("q0", "q1", "q2", "q3")

fixef_lst <- list()

for(i in 1:length(fits_phy)) {
  
  fit_phy <- fits_phy[i]
  
  load(paste0(dirPhy, fit_phy))
  
  fixef <- getFixef(fits = res$fits, robust = TRUE, 
                    estimates = res$R2_robust,
                    dimension = "phylogeny", Q = Qs[i], 
                    level = "all") 
  fixef_lst[[i]] <- fixef
  
}


fixef_table_phy <- do.call(rbind, fixef_lst)

##### By habitat #####
list.files("Results/Regressions/phylo-spec/Habitat/Scale/")

dirPhy_habitat <- "Results/Regressions/phylo-spec/Habitat/Scale/"

prefix <- "reg_Phylo_Spec_DIS_"

Qs <- c("q0", "q1", "q2", "q3")

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

fixef_q0_lst <- list()
fixef_q1_lst <- list()
fixef_q2_lst <- list()
fixef_q3_lst <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirPhy_habitat, prefix, Qs[1], "_", habitat[i], ".RData")) 
  
  fixefq0 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust,
                      dimension = "phylogeny", Q = Qs[1], 
                         level = habitat[i]) 
  fixef_q0_lst[[i]] <- fixefq0 
  
  # q1
  load(paste0(dirPhy_habitat, prefix, Qs[2], "_", habitat[i], ".RData")) 
  
  fixefq1 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust,
                      dimension = "phylogeny", Q = Qs[2], 
                         level = habitat[i]) 
  fixef_q1_lst[[i]] <- fixefq1
  
  # q2
  load(paste0(dirPhy_habitat, prefix, Qs[3], "_", habitat[i], ".RData")) 
  
  fixefq2 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust,
                      dimension = "phylogeny", Q = Qs[3], 
                         level = habitat[i]) 
  fixef_q2_lst[[i]] <- fixefq2
  
  # q3
  load(paste0(dirPhy_habitat, prefix, Qs[4], "_", habitat[i], ".RData")) 
  
  fixefq3 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust,
                      dimension = "phylogeny", Q = Qs[4], 
                         level = habitat[i]) 
  fixef_q3_lst[[i]] <- fixefq3
  
}

fixef_q0_table <- do.call(rbind, fixef_q0_lst)
fixef_q1_table <- do.call(rbind, fixef_q1_lst)
fixef_q2_table <- do.call(rbind, fixef_q2_lst)
fixef_q3_table <- do.call(rbind, fixef_q3_lst)

##### Combine all tables #####
library(tidyverse)

fixef_table_phylo <- rbind(fixef_table_phy, fixef_q0_table, fixef_q1_table, 
                        fixef_q2_table, fixef_q3_table)

write.csv(fixef_table_phylo, 
          file = "Results/Regressions/phylo-spec/Estimations_FIXEF_phylo.csv")

############### ------------ GET FIXEF TRAITS  ------------- ###############

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirTrt <- "Results/Regressions/trait-spec/SAM/Scale/"

fits_trt <- list.files(dirTrt)

Qs <- c("q0", "q1", "q2", "q3")

fixef_lst_dis <- list()
fixef_lst_met <- list()

for(i in 1:length(Qs)) {
  # Distance based metrics
  print(fits_trt[i])
  
  fit_trt_dis <- fits_trt[i]
  
  load(paste0(dirTrt, fit_trt_dis))
  
  fixefdis <- getFixef(fits = res$fits, robust = TRUE, 
                       estimates = res$R2_robust, 
                       dimension = "trait", Q = Qs[i], 
                       level = "all") 
  fixef_lst_dis[[i]] <- fixefdis
  
  # Classic metrics
  print(fits_trt[i + 4])
  
  fit_trt_met <- fits_trt[i + 4]
  
  load(paste0(dirTrt, fit_trt_met))
  
  fixefmet <- getFixef(fits = res$fits, robust = TRUE, 
                       estimates = res$R2_robust, 
                       dimension = "trait", Q = Qs[i], 
                       level = "all") 
  fixef_lst_met[[i]] <- fixefmet
  
}

fixef_table_trt_dis <- do.call(rbind, fixef_lst_dis)
fixef_table_trt_met <- do.call(rbind, fixef_lst_met)

##### By habitat #####
list.files("Results/Regressions/trait-spec/Habitat/DIS/")

dirTrt_habitat_dis <- "Results/Regressions/trait-spec/Habitat/DIS/"

prefixDIS <- "reg_Trait_Spec_DIS_SAM_scaled_"

Qs <- c("q0", "q1", "q2", "q3")

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

### Distance based metrics
fixef_q0_lst_dis <- list()
fixef_q1_lst_dis <- list()
fixef_q2_lst_dis <- list()
fixef_q3_lst_dis <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[1], "_", habitat[i], ".RData")) 
  
  fixefq0 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust, 
                      dimension = "trait", Q = Qs[1], 
                         level = habitat[i]) 
  fixef_q0_lst_dis[[i]] <- fixefq0 
  
  # q1
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[2], "_", habitat[i], ".RData")) 
  
  fixefq1 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust, 
                      dimension = "trait", Q = Qs[2], 
                         level = habitat[i]) 
  fixef_q1_lst_dis[[i]] <- fixefq1
  
  # q2
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[3], "_", habitat[i], ".RData")) 
  
  fixefq2 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust, 
                      dimension = "trait", Q = Qs[3], 
                         level = habitat[i]) 
  fixef_q2_lst_dis[[i]] <- fixefq2
  
  # q3
  load(paste0(dirTrt_habitat_dis, prefixDIS, Qs[4], "_", habitat[i], ".RData")) 
  
  fixefq3 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust, 
                      dimension = "trait", Q = Qs[4], 
                         level = habitat[i]) 
  fixef_q3_lst_dis[[i]] <- fixefq3
  
}

fixef_q0_table_dis <- do.call(rbind, fixef_q0_lst_dis)
fixef_q1_table_dis <- do.call(rbind, fixef_q1_lst_dis)
fixef_q2_table_dis <- do.call(rbind, fixef_q2_lst_dis)
fixef_q3_table_dis <- do.call(rbind, fixef_q3_lst_dis)

#### Classic metrics  
list.files("Results/Regressions/trait-spec/Habitat/MET/")

dirTrt_habitat_met <- "Results/Regressions/trait-spec/Habitat/MET/"

prefixMET <- "reg_Trait_Spec_MET_SAM_scaled_"

fixef_q0_lst_met <- list()
fixef_q1_lst_met <- list()
fixef_q2_lst_met <- list()
fixef_q3_lst_met <- list()

for(i in 1:length(habitat)) { 
  
  print(habitat[i])
  
  # q0
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[1], "_", habitat[i], ".RData")) 
  
  fixefq0 <- getFixef(fits = res$fits, robust = TRUE, 
                       estimates = res$R2_robust, 
                       dimension = "trait", Q = Qs[1], 
                         level = habitat[i]) 
  fixef_q0_lst_met[[i]] <- fixefq0 
  
  # q1
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[2], "_", habitat[i], ".RData")) 
  
  fixefq1 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust, 
                      dimension = "trait", Q = Qs[2], 
                         level = habitat[i]) 
  fixef_q1_lst_met[[i]] <- fixefq1
  
  # q2
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[3], "_", habitat[i], ".RData")) 
  
  fixefq2 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust, 
                      dimension = "trait", Q = Qs[3], 
                         level = habitat[i]) 
  fixef_q2_lst_met[[i]] <- fixefq2
  
  # q3
  load(paste0(dirTrt_habitat_met, prefixMET, Qs[4], "_", habitat[i], ".RData")) 
  
  fixefq3 <- getFixef(fits = res$fits, robust = TRUE, 
                      estimates = res$R2_robust, 
                      dimension = "trait", Q = Qs[4], 
                         level = habitat[i]) 
  fixef_q3_lst_met[[i]] <- fixefq3
  
}

fixef_q0_table_met <- do.call(rbind, fixef_q0_lst_met)
fixef_q1_table_met <- do.call(rbind, fixef_q1_lst_met)
fixef_q2_table_met <- do.call(rbind, fixef_q2_lst_met)
fixef_q3_table_met <- do.call(rbind, fixef_q3_lst_met)

##### Combine all tables #####
library(tidyverse)

fixef_table_trait_dis <- rbind(fixef_table_trt_dis, fixef_q0_table_dis, 
                            fixef_q1_table_dis, fixef_q2_table_dis, 
                            fixef_q3_table_dis)

fixef_table_trait_met <- rbind(fixef_table_trt_met, fixef_q0_table_met, 
                            fixef_q1_table_met, fixef_q2_table_met, 
                            fixef_q3_table_met)

write.csv(fixef_table_trait_dis, 
          file = "Results/Regressions/trait-spec/Estimations_FIXEF_trait_dis.csv")

write.csv(fixef_table_trait_met, 
          file = "Results/Regressions/trait-spec/Estimations_FIXEF_trait_met.csv")

############### ------------ GET FIXEF TAXONOMY ------------- ###############

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirTax <- "Results/Regressions/taxo-spec/"

fits_tax <- list.files(dirTax)

load(paste0(dirTax, "reg_Taxo_Spec_Alpha_THRESH.RData"))
  
fixeftax <- getFixef(fits = res$fits, robust = TRUE, 
                     estimates = res$R2_robust, 
                     dimension = "taxonomy", Q = "Q13", 
                     level = "all") 

##### By habitat #####

dirTaxH <- "Results/Regressions/taxo-spec/Habitat/"

fits_tax_habitat <- list.files(dirTaxH)

habitat <- c("deciduousForest", "evergreenForest", "grasslandHerbaceous", 
             "mixedForest", "shrubScrub")

fixef_lst_tax <- list()

for(i in 1:length(habitat)) {
  
  fit <- fits_tax_habitat[6:10][i] 
  
  load(paste0(dirTaxH, fit))
  
  fixeftaxH <- getFixef(fits = res$fits, robust = TRUE, 
                        estimates = res$R2_robust, 
                        dimension = "taxonomy", Q = "Q13", 
                        level = habitat[i]) 
  
  fixef_lst_tax[[i]] <- fixeftaxH 
  
}

fixef_table_habitat <- do.call(rbind, fixef_lst_tax)

fixef_table_taxonomy <- rbind(fixeftax, fixef_table_habitat)

write.csv(fixef_table_taxonomy, 
          file = "Results/Regressions/taxo-spec/Estimations_FIXEF_taxonomy.csv")
