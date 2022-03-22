
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
                             nDraws = 1000, 
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
                                 nDraws = 1000, 
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
                              nDraws = 1000, 
                              dimension = "trait", 
                              Q = Qs[1], 
                              level = habitat[i]) 
  beta_q0_lst_dis[[i]] <- betaq0 
  
  # q1
  load(paste0(dirTrt_habitat_dis, prefix_dis, Qs[2], "_", habitat[i], ".RData")) 
  
  betaq1 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
                              dimension = "trait",
                              Q = Qs[2], 
                              level = habitat[i]) 
  beta_q1_lst_dis[[i]] <- betaq1
  
  # q2
  load(paste0(dirTrt_habitat_dis, prefix_dis, Qs[3], "_", habitat[i], ".RData")) 
  
  betaq2 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
                              dimension = "trait",
                              Q = Qs[3], 
                              level = habitat[i]) 
  beta_q2_lst_dis[[i]] <- betaq2
  
  # q3
  load(paste0(dirTrt_habitat_dis, prefix_dis, Qs[4], "_", habitat[i], ".RData")) 
  
  betaq3 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
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
                              nDraws = 1000, 
                              dimension = "trait", 
                              Q = Qs[1], 
                              level = habitat[i]) 
  beta_q0_lst_met[[i]] <- betaq0 
  
  # q1
  load(paste0(dirTrt_habitat_met, prefix_met, Qs[2], "_", habitat[i], ".RData")) 
  
  betaq1 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
                              dimension = "trait",
                              Q = Qs[2], 
                              level = habitat[i]) 
  beta_q1_lst_met[[i]] <- betaq1
  
  # q2
  load(paste0(dirTrt_habitat_met, prefix_met, Qs[3], "_", habitat[i], ".RData")) 
  
  betaq2 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
                              dimension = "trait",
                              Q = Qs[3], 
                              level = habitat[i]) 
  beta_q2_lst_met[[i]] <- betaq2
  
  # q3
  load(paste0(dirTrt_habitat_met, prefix_met, Qs[4], "_", habitat[i], ".RData")) 
  
  betaq3 <- getFixefPosterior(fits = res$fits, 
                              estimates = res$R2_robust, 
                              nDraws = 1000, 
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
