### Libraries
library(tidyverse)
library(tidymodels)

### Auxiliary function
source("R/NEON_diversity/R/Functions/demonPRED.R")

########## --------- Predictions phylogenetic dimension under distance matrices ---------- ##########

##### Load and prepare data ######
## Q0
phylo_spec_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q0.csv") 

phylo_spec_NEON_table_q0

## Q1
phylo_spec_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q1.csv")

phylo_spec_NEON_table_q1

## Q2
phylo_spec_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q2.csv")

phylo_spec_NEON_table_q2

### data partition

### Initial splits
## Q0
set.seed(12345)

data_split_q0 <- phylo_spec_NEON_table_q0 %>% 
  initial_split(prop = 0.50) # proportion of data to be retained for modeling 
## Training data
data_train_q0 <- training(data_split_q0)
## Testing data
data_test_q0 <- testing(data_split_q0)

## Q1
set.seed(12345)

data_split_q1 <- phylo_spec_NEON_table_q1  %>% 
  initial_split(prop = 0.50) # proportion of data to be retained for modeling 
## Training data
data_train_q1 <- training(data_split_q1)
## Testing data
data_test_q1 <- testing(data_split_q1)

## Q2
set.seed(12345)

data_split_q2 <- phylo_spec_NEON_table_q2 %>% 
  initial_split(prop = 0.50) # proportion of data to be retained for modeling 
## Training data
data_train_q2 <- training(data_split_q2)
## Testing data
data_test_q2 <- testing(data_split_q2)

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

### Run PD
pred_PD_phylo <- godPrediction(trainData = data_train_q0,  # data to train the model 
                               testData = data_test_q0, # data to test the model
                               variable = "PD_phylo", 
                               covariable = "SD_spec", 
                               dimension = "phylogeny", 
                               Q = Qs[1],
                               nCores = nCores, 
                               nIters = nIters, 
                               burnin = burnin, 
                               nChains = nChains, 
                               control = control, 
                               engine = engine
) 

### Run PDz
pred_PDz_phylo <- godPrediction(trainData = data_train_q0,  # data to train the model 
                                testData = data_test_q0, # data to test the model
                                variable = "PDz_phylo", 
                                covariable = "SD_spec", 
                                dimension = "phylogeny", 
                                Q = Qs[1],
                                nCores = nCores, 
                                nIters = nIters, 
                                burnin = burnin, 
                                nChains = nChains, 
                                control = control, 
                                engine = engine
) 

### Run MPD
pred_MPD_phylo <- godPrediction(trainData = data_train_q0,  # data to train the model 
                                testData = data_test_q0, # data to test the model
                                variable = "MPD_phylo", 
                                covariable = "MSD_spec", 
                                dimension = "phylogeny", 
                                Q = Qs[1],
                                nCores = nCores, 
                                nIters = nIters, 
                                burnin = burnin, 
                                nChains = nChains, 
                                control = control, 
                                engine = engine
) 

### Run MPDz
pred_MPDz_phylo <- godPrediction(trainData = data_train_q0,  # data to train the model 
                                 testData = data_test_q0, # data to test the model
                                 variable = "MPDz_phylo", 
                                 covariable = "MSD_spec", 
                                 dimension = "phylogeny", 
                                 Q = Qs[1],
                                 nCores = nCores, 
                                 nIters = nIters, 
                                 burnin = burnin, 
                                 nChains = nChains, 
                                 control = control, 
                                 engine = engine
) 

### Run q0_DPM
pred_qDTM_phylo_DIS_q0 <- godPrediction(trainData = data_train_q0,  # data to train the model 
                                        testData = data_test_q0, # data to test the model
                                        variable = "qDTM_phylo", 
                                        covariable = "qDTM_spec", 
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
pred_qDTM_phylo_DIS_q1 <- godPrediction(trainData = data_train_q1,  # data to train the model 
                                        testData = data_test_q1, # data to test the model
                                        variable = "qDTM_phylo", 
                                        covariable = "qDTM_spec", 
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
pred_qDTM_phylo_DIS_q2 <- godPrediction(trainData = data_train_q2,  # data to train the model 
                                        testData = data_test_q2, # data to test the model
                                        variable = "qDTM_phylo", 
                                        covariable = "qDTM_spec", 
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
pred_phylo_DIS <- list(pred_PD_phylo, 
                       pred_PDz_phylo, 
                       pred_MPD_phylo, 
                       pred_MPDz_phylo, 
                       pred_qDTM_phylo_DIS_q0, 
                       pred_qDTM_phylo_DIS_q1, 
                       pred_qDTM_phylo_DIS_q2)

### Save results
save(pred_phylo_DIS,  
     file = "output/Predictions/Reanalyses/predictions_phylo_DIS.RData")

########## --------- Predictions trait dimension under distance matrices ---------- ##########

##### Load and prepare data ######
## Q0
trait_spec_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q0.csv") 

trait_spec_NEON_table_q0

## Q1
trait_spec_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q1.csv")

trait_spec_NEON_table_q1

## Q2
trait_spec_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q2.csv")

trait_spec_NEON_table_q2

### data partition

### Initial splits
## Q0
set.seed(12345)

data_split_q0 <- trait_spec_NEON_table_q0 %>% 
  initial_split(prop = 0.50) # proportion of data to be retained for modeling 
## Training data
data_train_q0 <- training(data_split_q0)
## Testing data
data_test_q0 <- testing(data_split_q0)

## Q1
set.seed(12345)

data_split_q1 <- trait_spec_NEON_table_q1  %>% 
  initial_split(prop = 0.50) # proportion of data to be retained for modeling 
## Training data
data_train_q1 <- training(data_split_q1)
## Testing data
data_test_q1 <- testing(data_split_q1)

## Q2
set.seed(12345)

data_split_q2 <- trait_spec_NEON_table_q2 %>% 
  initial_split(prop = 0.50) # proportion of data to be retained for modeling 
## Training data
data_train_q2 <- training(data_split_q2)
## Testing data
data_test_q2 <- testing(data_split_q2)

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

### Run SR
pred_SR_trait <- godPrediction(trainData = data_train_q0,  # data to train the model 
                               testData = data_test_q0, # data to test the model
                               variable = "SR_trait", 
                               covariable = "MSD_spec", 
                               dimension = "trait", 
                               Q = Qs[1],
                               nCores = nCores, 
                               nIters = nIters, 
                               burnin = burnin, 
                               nChains = nChains, 
                               control = control, 
                               engine = engine
) 

### Run MTD
pred_MTD_trait <- godPrediction(trainData = data_train_q0,  # data to train the model 
                                testData = data_test_q0, # data to test the model
                                variable = "MFD_trait", 
                                covariable = "MSD_spec", 
                                dimension = "trait", 
                                Q = Qs[1],
                                nCores = nCores, 
                                nIters = nIters, 
                                burnin = burnin, 
                                nChains = nChains, 
                                control = control, 
                                engine = engine
) 

### Run MTDz
pred_MTDz_trait <- godPrediction(trainData = data_train_q0,  # data to train the model 
                                 testData = data_test_q0, # data to test the model
                                 variable = "MFDz_trait", 
                                 covariable = "MSD_spec", 
                                 dimension = "trait", 
                                 Q = Qs[1],
                                 nCores = nCores, 
                                 nIters = nIters, 
                                 burnin = burnin, 
                                 nChains = nChains, 
                                 control = control, 
                                 engine = engine
) 

### Run q0_DTM
pred_qDTM_trait_DIS_q0 <- godPrediction(trainData = data_train_q0,  # data to train the model 
                                        testData = data_test_q0, # data to test the model
                                        variable = "qDTM_trait", 
                                        covariable = "qDTM_spec", 
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
pred_qDTM_trait_DIS_q1 <- godPrediction(trainData = data_train_q1,  # data to train the model 
                                        testData = data_test_q1, # data to test the model
                                        variable = "qDTM_trait", 
                                        covariable = "qDTM_spec", 
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
pred_qDTM_trait_DIS_q2 <- godPrediction(trainData = data_train_q2,  # data to train the model 
                                        testData = data_test_q2, # data to test the model
                                        variable = "qDTM_trait", 
                                        covariable = "qDTM_spec", 
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
pred_trait_DIS <- list(pred_SR_trait, 
                       pred_MTD_trait, 
                       pred_MTDz_trait, 
                       pred_qDTM_trait_DIS_q0, 
                       pred_qDTM_trait_DIS_q1, 
                       pred_qDTM_trait_DIS_q2)

### Save results
save(pred_trait_DIS,  
     file = "output/Predictions/Reanalyses/predictions_trait_DIS.RData")
