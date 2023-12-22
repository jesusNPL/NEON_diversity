### Libraries
library(tidyverse)
library(tidymodels)

### Auxiliary function
source("R/NEON_diversity/R/Functions/demonPRED.R")

########## --------- Predictions phylogenetic dimension under SS ---------- ##########

##### Load and prepare data ######
phylo_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_SS_NEON.csv") 

phylo_spec_SS_NEON

### data partition
set.seed(12345)

### Initial splits
data_split <- phylo_spec_SS_NEON %>% 
  initial_split(prop = 0.50) # proportion of data to be retained for modeling 

## Training data
data_train <- training(data_split)
## Testing data
data_test <- testing(data_split)

### Inits 
# Diversity orders
Qs <- c("q0", "q1", "q2") 
# Number of cores 
nCores <- 20 
# Number of generations
nIters <- 5000 
# Burnin 20%
burnin <- nIters/5 
# number of NUTS
nChains <- 4 
# Additional control
control <- list(adapt_delta = 0.99) 
# backend 
engine <- "cmdstanr"

##### Run model predictions #####

### Run q0_DPM
pred_qDTM_phylo_q0 <- godPrediction(trainData = data_train,  # data to train the model 
                                    testData = data_test, # data to test the model
                                    variable = "qDTM_phylo_q0", 
                                    covariable = "SRich_spec", 
                                    dimension = "phylogeny", 
                                    Q = Qs[1],
                                    nCores = nCores, 
                                    nIters = nIters, 
                                    burnin = burnin, 
                                    nChains = nChains, 
                                    control = control, 
                                    engine = engine
) 

### Run q1_DPM
pred_qDTM_phylo_q1 <- godPrediction(trainData = data_train,  # data to train the model 
                                    testData = data_test, # data to test the model
                                    variable = "qDTM_phylo_q1", 
                                    covariable = "Shannon_spec", 
                                    dimension = "phylogeny", 
                                    Q = Qs[2],
                                    nCores = nCores, 
                                    nIters = nIters, 
                                    burnin = burnin, 
                                    nChains = nChains, 
                                    control = control, 
                                    engine = engine
) 

### Run q2_DPM
pred_qDTM_phylo_q2 <- godPrediction(trainData = data_train,  # data to train the model 
                                    testData = data_test, # data to test the model
                                    variable = "qDTM_phylo_q2", 
                                    covariable = "Simpson_spec", 
                                    dimension = "phylogeny", 
                                    Q = Qs[3],
                                    nCores = nCores, 
                                    nIters = nIters, 
                                    burnin = burnin, 
                                    nChains = nChains, 
                                    control = control, 
                                    engine = engine
) 

## List of results
pred_phylo_SS <- list(pred_qDTM_phylo_q0, 
                      pred_qDTM_phylo_q1, 
                      pred_qDTM_phylo_q2)

### Save results
save(pred_phylo_SS,  
     file = "output/Predictions/Reanalyses/predictions_phylo_SS.RData")

########## --------- Predictions trait dimension under SS ---------- ##########

##### Load and prepare data ######
trait_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/trait_spec_SS_NEON.csv") 

trait_spec_SS_NEON

### data partition
set.seed(12345)

### Initial splits
data_split <- trait_spec_SS_NEON %>% 
  initial_split(prop = 0.50) # proportion of data to be retained for modeling 

## Training data
data_train <- training(data_split)
## Testing data
data_test <- testing(data_split)

### Inits
# Diversity orders
Qs <- c("q0", "q1", "q2") 
# Number of cores 
nCores <- 20 
# Number of generations
nIters <- 5000 
# Burnin 20%
burnin <- nIters/5 
# number of NUTS
nChains <- 4 
# Additional control
control <- list(adapt_delta = 0.99) 
# backend 
engine <- "cmdstanr"

##### Run model predictions #####

### Run q0_DTM
pred_qDTM_trait_q0 <- godPrediction(trainData = data_train,  # data to train the model 
                                    testData = data_test, # data to test the model
                                    variable = "qDTM_trait_q0", 
                                    covariable = "SRich_spec", 
                                    dimension = "trait", 
                                    Q = Qs[1],
                                    nCores = nCores, 
                                    nIters = nIters, 
                                    burnin = burnin, 
                                    nChains = nChains, 
                                    control = control, 
                                    engine = engine
) 

### Run q1_DTM
pred_qDTM_trait_q1 <- godPrediction(trainData = data_train,  # data to train the model 
                                    testData = data_test, # data to test the model
                                    variable = "qDTM_trait_q1", 
                                    covariable = "Shannon_spec", 
                                    dimension = "trait", 
                                    Q = Qs[2],
                                    nCores = nCores, 
                                    nIters = nIters, 
                                    burnin = burnin, 
                                    nChains = nChains, 
                                    control = control, 
                                    engine = engine
) 

### Run q2_DTM
pred_qDTM_trait_q2 <- godPrediction(trainData = data_train,  # data to train the model 
                                    testData = data_test, # data to test the model
                                    variable = "qDTM_trait_q2", 
                                    covariable = "Simpson_spec", 
                                    dimension = "trait", 
                                    Q = Qs[3],
                                    nCores = nCores, 
                                    nIters = nIters, 
                                    burnin = burnin, 
                                    nChains = nChains, 
                                    control = control, 
                                    engine = engine
) 

## List of results
pred_trait_SS <- list(pred_qDTM_trait_q0, 
                      pred_qDTM_trait_q1, 
                      pred_qDTM_trait_q2)

### Save results
save(pred_trait_SS,  
     file = "output/Predictions/Reanalyses/predictions_trait_SS.RData")

########## --------- Predictions taxonomic dimension under SS ---------- ##########

##### Load and prepare data ######
taxo_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/taxo_spec_SS_NEON.csv") %>% 
  mutate(siteID = Site) %>% 
  select(!Site) %>% 
  select(siteID, everything())

taxo_spec_SS_NEON

### data partition
set.seed(12345)

### Initial splits
data_split <- taxo_spec_SS_NEON %>% 
  initial_split(prop = 0.50) # proportion of data to be retained for modeling 

## Training data
data_train <- training(data_split)
## Testing data
data_test <- testing(data_split)

### Inits
# Diversity orders
Qs <- c("q0", "q1", "q2") 
# Number of cores 
nCores <- 20 
# Number of generations
nIters <- 5000 
# Burnin 20%
burnin <- nIters/5 
# number of NUTS
nChains <- 4 
# Additional control
control <- list(adapt_delta = 0.99) 
# backend 
engine <- "cmdstanr"

##### Run model predictions #####

### Run q0_S
pred_S_taxo_q0 <- godPrediction(trainData = data_train,  # data to train the model 
                                testData = data_test, # data to test the model
                                variable = "SRich_ground", 
                                covariable = "SRich_spec", 
                                dimension = "taxonomy", 
                                Q = Qs[1],
                                nCores = nCores, 
                                nIters = nIters, 
                                burnin = burnin, 
                                nChains = nChains, 
                                control = control, 
                                engine = engine
) 

### Run q1_H
pred_H_taxo_q1 <- godPrediction(trainData = data_train,  # data to train the model 
                                testData = data_test, # data to test the model
                                variable = "Shannon_ground", 
                                covariable = "Shannon_spec", 
                                dimension = "taxonomy", 
                                Q = Qs[2],
                                nCores = nCores, 
                                nIters = nIters, 
                                burnin = burnin, 
                                nChains = nChains, 
                                control = control, 
                                engine = engine
) 

### Run q2_D
pred_D_taxo_q2 <- godPrediction(trainData = data_train,  # data to train the model 
                                testData = data_test, # data to test the model
                                variable = "Simpson_ground", 
                                covariable = "Simpson_spec", 
                                dimension = "taxonomy", 
                                Q = Qs[3],
                                nCores = nCores, 
                                nIters = nIters, 
                                burnin = burnin, 
                                nChains = nChains, 
                                control = control, 
                                engine = engine
) 

## List of results
pred_taxo_SS <- list(pred_S_taxo_q0, 
                      pred_H_taxo_q1, 
                      pred_D_taxo_q2)

### Save results
save(pred_taxo_SS,  
     file = "output/Predictions/Reanalyses/predictions_taxo_SS.RData")

