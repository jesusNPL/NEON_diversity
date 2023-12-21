library(tidyverse)
### Auxiliary function
source("R/NEON_diversity/R/Functions/GOD_Bayes_taxo_spec.R")

########## ---------- Taxonomic dimension ---------- ##########

##### Load and prepare data for BLMM - taxonomic dimension ######
taxo_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/taxo_spec_SS_NEON.csv") 

taxo_spec_SS_NEON <- taxo_spec_SS_NEON %>% 
  select(siteID = Site, # group level variable 
         plotID, 
         SRich_ground, Shannon_ground, Simpson_ground, # ground based metrics
         SRich_spec, Shannon_spec, Simpson_spec) # spectral based metrics

taxo_spec_SS_NEON

##### Bayesian multilevel model between taxonomic dimension and spectral dimension based on SS #####

### Inits
# Number of metrics
nMetrics <- 3
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
# Path to save the results
saveFilePath <- "Results/Regressions/Reanalyses/taxo_SS/BLM_taxo_SS_alpha_scale.RData"

### Run BLMM
god_BayReg_alpha_taxo(resMetrics = taxo_spec_SS_NEON, 
                      nMetrics = nMetrics, 
                      scale = TRUE, 
                      nChains = nChains, 
                      nIters = nIters, 
                      burnin = burnin, 
                      nCores = nCores, 
                      engine = engine, 
                      pathSave = saveFilePath
                      )

### Run BLMM not scaling SS

# Path to save the results
saveFilePath <- "Results/Regressions/Reanalyses/taxo_SS/BLM_taxo_SS_alpha.RData"

god_BayReg_alpha_taxo(resMetrics = taxo_spec_SS_NEON, 
                      nMetrics = nMetrics, 
                      scale = FALSE, 
                      nChains = nChains, 
                      nIters = nIters, 
                      burnin = burnin, 
                      nCores = nCores, 
                      engine = engine, 
                      pathSave = saveFilePath
)

########## ---------- Phylogenetic dimension ---------- ##########

### Auxiliary function
source("R/NEON_diversity/R/Functions/GOD_Bayes_taxo_spec.R")

##### Load and prepare data for BLMM - Phylogenetic dimension ######
phylo_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_SS_NEON.csv") 

phylo_spec_SS_NEON

##### Bayesian multilevel model between phylogenetic dimension and spectral dimension based on SS #####

### Inits
# Number of metrics
nMetrics <- 3
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
# Path to save the results
saveFilePath <- "Results/Regressions/Reanalyses/taxo_SS/BLM_phylo_SS_alpha_scale.RData"

### Run BLMM
god_BayReg_alpha_taxo(resMetrics = phylo_spec_SS_NEON, 
                      nMetrics = nMetrics, 
                      scale = TRUE, 
                      nChains = nChains, 
                      nIters = nIters, 
                      burnin = burnin, 
                      nCores = nCores, 
                      engine = engine, 
                      pathSave = saveFilePath
)

### Run BLMM not scaling SS
# Path to save the results
saveFilePath <- "Results/Regressions/Reanalyses/taxo_SS/BLM_phylo_SS_alpha.RData"

god_BayReg_alpha_taxo(resMetrics = phylo_spec_SS_NEON, 
                      nMetrics = nMetrics, 
                      scale = FALSE, 
                      nChains = nChains, 
                      nIters = nIters, 
                      burnin = burnin, 
                      nCores = nCores, 
                      engine = engine, 
                      pathSave = saveFilePath
)

########## ---------- Trait dimension ---------- ##########

### Auxiliary function
source("R/NEON_diversity/R/Functions/GOD_Bayes_taxo_spec.R")

##### Load and prepare data for BLMM - Phylogenetic dimension ######
trait_spec_SS_NEON <- read_csv("Results/RegDATA/Reanalyses/trait_spec_SS_NEON.csv") 

trait_spec_SS_NEON

##### Bayesian multilevel model between traitgenetic dimension and spectral dimension based on SS #####

### Inits
# Number of metrics
nMetrics <- 3
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
# Path to save the results
saveFilePath <- "Results/Regressions/Reanalyses/taxo_SS/BLM_trait_SS_alpha_scale.RData"

### Run BLMM
god_BayReg_alpha_taxo(resMetrics = trait_spec_SS_NEON, 
                      nMetrics = nMetrics, 
                      scale = TRUE, 
                      nChains = nChains, 
                      nIters = nIters, 
                      burnin = burnin, 
                      nCores = nCores, 
                      engine = engine, 
                      pathSave = saveFilePath
)

### Run BLMM not scaling SS
# Path to save the results
saveFilePath <- "Results/Regressions/Reanalyses/taxo_SS/BLM_trait_SS_alpha.RData"

god_BayReg_alpha_taxo(resMetrics = trait_spec_SS_NEON, 
                      nMetrics = nMetrics, 
                      scale = FALSE, 
                      nChains = nChains, 
                      nIters = nIters, 
                      burnin = burnin, 
                      nCores = nCores, 
                      engine = engine, 
                      pathSave = saveFilePath
)

