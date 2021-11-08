###### Part 1 - run taxonomic dimension #####

library(dplyr)
library(tidyr)
library(picante)
# remotes::install_github("cran/vegetarian")
library(vegetarian)
library(vegan)
source("R/NEON_diversity/R/Functions/SS_diversity.R")
source("R/NEON_diversity/R/Functions/auxiliar.R")

load("DATA/RData/HARV_SS.RData") # Spectral species
load("DATA/Harvard_NEON/Data/harv_matchData.RData") # HARV data
rm(matchPhylo, matchTrait)

#### Prepare data #####

# Spectral species with threshold
specSPP_thresh <- paste0("SSpp_", ss_threshold$SS)
ss_threshold$specSpp <- specSPP_thresh

ss_threshold <- ss_threshold[, c(1:4, 431, 432, 5:430)]
ss_threshold[1:5, 1:7]

# Spectral species without threshold
specSPP_no_thresh <- paste0("SSpp_", ss_NO_threshold$SS)
ss_NO_threshold$specSpp <- specSPP_no_thresh

ss_NO_threshold <- ss_NO_threshold[, c(1:4, 431, 432, 5:430)]
ss_NO_threshold[1:5, 1:7]

### Make community data matrices
harv_commSPEC_thresh <- make_SSppComm(spectra = ss_threshold)

harv_commSPEC_no_thresh <- make_SSppComm(spectra = ss_NO_threshold)

### Clean vegetation community data

harv_comm <- harv_comm[harv_comm$plotID_SPEC %in% rownames(harv_commSPEC_no_thresh), ]
rownames(harv_comm) <- harv_comm$plotID_SPEC
harv_comm <- harv_comm[match(rownames(harv_commSPEC_no_thresh), harv_comm$plotID_SPEC), ]

harv_comm[1:5, 1:7]

harv_comm_clean <- harv_comm_clean[harv_comm_clean$plotID_SPEC %in% rownames(harv_commSPEC_no_thresh), ]
rownames(harv_comm_clean) <- harv_comm_clean$plotID_SPEC
# Order plots according plots in spectra
harv_comm_clean <- harv_comm_clean[match(rownames(harv_commSPEC_no_thresh), harv_comm_clean$plotID_SPEC), ]

harv_comm_clean[1:5, 1:5]

##### Spectral Species diversity under different q #####
plotNames <- rownames(harv_commSPEC_no_thresh)
nPlots <- length(plotNames)

spec_HARV_div_thresh <- demon_SSppDIV(
  comm = harv_commSPEC_thresh, q1 = 0, q2 = 1, q3 = 2,
  plotNames = plotNames, nPlots = nPlots
)

spec_HARV_div_no_thresh <- demon_SSppDIV(
  comm = harv_commSPEC_no_thresh, q1 = 0, q2 = 1, q3 = 2,
  plotNames = plotNames, nPlots = nPlots
)

##### Vegetation diversity under different q #####

veg_HARV_div <- demon_SSppDIV(
  comm = harv_comm_clean[, 4:121], q1 = 0, q2 = 1, q3 = 2,
  plotNames = plotNames, nPlots = nPlots
)

##### Combine information  of the two dimensions #####
Div_tax_alpha_thresh <- makeTable_taxonomy(
  div1 = veg_HARV_div$commDiv,
  div2 = spec_HARV_div_thresh$commDiv
)

Div_tax_alpha_thresh$divOBS$herv_cov <- harv_comm$herb_cov
Div_tax_alpha_thresh$divOBS$tree_cov <- harv_comm$tree_cov

Div_tax_alpha_no_thresh <- makeTable_taxonomy(
  div1 = veg_HARV_div$commDiv,
  div2 = spec_HARV_div_no_thresh$commDiv
)

Div_tax_alpha_no_thresh$divOBS$herv_cov <- harv_comm$herb_cov
Div_tax_alpha_no_thresh$divOBS$tree_cov <- harv_comm$tree_cov

### Save tables with alpha calculations
saveRDS(Div_tax_alpha_no_thresh, file = "Results/HARV/taxo_alpha_no_thresh.rds")
saveRDS(Div_tax_alpha_thresh, file = "Results/HARV/taxo_alpha_thresh.rds")

save(Div_tax_alpha_thresh, Div_tax_alpha_no_thresh,
  file = "Results/HARV/div_tax_harv.RData"
)

###### Part 2 - run phylogenetic dimension #####

##### Run alpha diversity Phylogenetic and functional #####

library(dplyr)
library(tidyr)
library(picante)
library(fundiversity)
library(phytools)
source("R/NEON_diversity/R/Functions/AlphaSpecDIV.R")

load("DATA/RData/HARV_SS.RData") # Spectral species
rm(ss_NO_threshold)

## Spectra data range from column 4 to 430
plotNames <- unique(ss_threshold$plotID)
nPlots <- length(plotNames)
Range <- 5:430 ## Spectra data range from column 5 to 430
QS <- c(0, 1, 2, 3) # different q values

spec_harv_dis <- list()

for (i in 1:length(QS)) {
  print(paste0("Starting calculations under ", QS[i], "..."))

  spec_harv_dis[[i]] <- demon_DivDistance(
    spectra = ss_threshold, plotNames = plotNames,
    nPlots = nPlots, specRange = Range, Q = QS[i],
    standardize = FALSE, gowDist = TRUE
  )
}

spec_harv_table <- do.call(rbind, spec_harv_dis)

### save results
save(spec_harv_dis, spec_harv_table, file = "Results/HARV/div_spec_harv.RData")

##### Run phylogenetic diversity #####

load("DATA/Harvard_NEON/Data/harv_matchData.RData") # HARV data
rm(harv_comm, harv_taxonomy, matchTrait)

harvNames <- unique(harv_comm_clean$plotID_SPEC)
QS <- c(0, 1, 2, 3) # different q values

phylo_harv_dis <- list()

for (j in 1:length(QS)) {
  print(paste0("Starting calculations under q = ", QS[j], "..."))

  phylo_harv_dis[[j]] <- demon_PhyloDistance(
    comm = matchPhylo$comm,
    phylo = matchPhylo$phy,
    abundance = TRUE, Q = QS[j],
    plotNames = harvNames
  )
}

phylo_harv_table <- do.call(rbind, phylo_harv_dis)

### save results
save(phylo_harv_dis, phylo_harv_table, file = "Results/HARV/div_phylo_harv.RData")

##### Combine phylogenetic and spectral distance results #####
Div_phylo_spec_alpha <- makeTable_phylogeny(
  div1 = phylo_harv_table,
  div2 = spec_harv_table
)

save(Div_phylo_spec_alpha, file = "Results/HARV/div_phylo_spec_harv.RData")

###### Part 3 - run trait dimension #####

##### Run functional diversity #####
library(dplyr)
library(tidyr)
library(picante)
library(fundiversity)
library(phytools)
library(gawdis)
source("R/NEON_diversity/R/Functions/AlphaSpecDIV.R")

load("DATA/Harvard_NEON/Data/harv_matchData.RData") # HARV data
rm(harv_comm, harv_taxonomy)

# trait_BEIN_NEON <- readRDS("DATA/Traits/Traits_BIEN/traits_BIEN_NEON_spp_genus_combined.rds")
# matchTrait <- trait_BEIN_NEON
# rownames(matchTrait) <- matchTrait$Taxa
# matchTrait <- match.phylo.data(phy = matchPhylo$phy, matchTrait)$data
# count no NAs
apply(matchTrait, MARGIN = 2, function(x) sum(!is.na(x)))

## Select traits with at least 50 observations

matchTrait_sel <- matchTrait %>%
  select(
    Taxa, mean_WPH, mean_MaxWPH, mean_LNCPLDryMass, mean_LArea,
    mean_LAreaPLDryMass, mean_SeedMass, mean_LDryMass
  ) %>%
  drop_na(mean_WPH)
# count no NAs
apply(matchTrait_sel, MARGIN = 2, function(x) sum(!is.na(x)))
# check if the values are numeric
apply(matchTrait_sel, MARGIN = 2, class)
# transform to numeric
matchTrait_sel[, 2:8] <- lapply(matchTrait_sel[, 2:8], function(x) as.numeric(as.character(x)))
rownames(matchTrait_sel) <- matchTrait_sel$Taxa

# Select only species present in the trait database
harv_comm <- t(harv_comm_clean[, 4:121])
harv_comm <- harv_comm[rownames(harv_comm) %in% rownames(matchTrait_sel), ]
harv_comm <- t(harv_comm)

### Run trait diversity
plotNames <- harv_comm_clean$plotID_SPEC
nPlots <- length(plotNames)

QS <- c(0, 1, 2, 3) # different q values

trait_harv_dis <- list()

for (i in 1:length(QS)) {
  print(paste0("Sarting calculations under Q = ", QS[i]))

  trait_harv_dis[[i]] <- demon_TraitDistance(
    comm = harv_comm, trait = matchTrait_sel[, 2:8], abundance = TRUE,
    Q = QS[i], plotNames = plotNames,
    nPlots = nPlots, gowDist = TRUE, gawDist = FALSE
  )
}

trait_harv_table <- do.call(rbind, trait_harv_dis)

### save results
save(trait_harv_dis, trait_harv_table, file = "Results/HARV/div_trait_harv.RData")

##### Run spectral diversity #####
load("DATA/RData/HARV_SS.RData") # Spectral species
rm(ss_NO_threshold)

## Spectra data range from column 4 to 430
plotNames <- unique(ss_threshold$plotID)
nPlots <- length(plotNames)
Range <- 5:430 ## Spectra data range from column 5 to 430
QS <- c(0, 1, 2, 3) # different q values

spec_harv_dis <- list()

for (i in 1:length(QS)) {
  print(paste0("Starting calculations under ", QS[i], "..."))

  spec_harv_dis[[i]] <- demon_DivDistance(
    spectra = ss_threshold, plotNames = plotNames,
    nPlots = nPlots, specRange = Range, Q = QS[i],
    standardize = FALSE, gowDist = TRUE
  )
}

spec_harv_table <- do.call(rbind, spec_harv_dis)

### save results
save(spec_harv_dis, spec_harv_table, file = "Results/HARV/div_spec_harv.RData")

##### Combine phylogenetic and spectral distance results #####
Div_trait_spec_alpha <- makeTable_traits(
  div1 = trait_harv_table,
  div2 = spec_harv_table
)

save(Div_trait_spec_alpha, file = "Results/HARV/div_trait_spec_harv.RData")

###### Part 4 - Bayesian correlations #####

##### Make Bayesian correlations taxonomy #####
library(tidyverse)

source("R/NEON_diversity/R/Functions/BayesianCorrelation.R")

load("Results/HARV/div_tax_spec_harv.RData")

corr_harv_tax_thresh <- demon_BayCorrTAX(
  matRES = Div_tax_alpha_thresh$divOBS,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_TAX_thresh_correlations.RData"
)

corr_harv_tax_no_thresh <- demon_BayCorrTAX(
  matRES = Div_tax_alpha_no_thresh$divOBS,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_TAX_no_thresh_correlations.RData"
)

##### Make Bayesian correlations phylogeny #####

load("Results/HARV/div_phylo_spec_harv.RData")

corr_harv_PHY_Q0 <- demon_BayCorrPHY(
  matRES = Div_phylo_spec_alpha,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_Phylo_Q0_correlations.RData",
  Q = 0
)

corr_harv_PHY_Q1 <- demon_BayCorrPHY(
  matRES = Div_phylo_spec_alpha,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_Phylo_Q1_correlations.RData",
  Q = 1
)

corr_harv_PHY_Q2 <- demon_BayCorrPHY(
  matRES = Div_phylo_spec_alpha,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_Phylo_Q2_correlations.RData",
  Q = 2
)

corr_harv_PHY_Q3 <- demon_BayCorrPHY(
  matRES = Div_phylo_spec_alpha,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_Phylo_Q3_correlations.RData",
  Q = 3
)

##### Make Bayesian correlations traits #####
source("R/NEON_diversity/R/Functions/BayesianCorrelation.R")

load("Results/HARV/div_trait_spec_harv.RData")

corr_harv_TRAIT_Q0 <- demon_BayCorrTRAIT(
  matRES = Div_trait_spec_alpha,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_Trait_Q0_correlations.RData",
  Q = 0
)

corr_harv_TRAIT_Q1 <- demon_BayCorrTRAIT(
  matRES = Div_trait_spec_alpha,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_Trait_Q1_correlations.RData",
  Q = 1
)

corr_harv_TRAIT_Q2 <- demon_BayCorrTRAIT(
  matRES = Div_trait_spec_alpha,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_Trait_Q2_correlations.RData",
  Q = 2
)

corr_harv_TRAIT_Q3 <- demon_BayCorrTRAIT(
  matRES = Div_trait_spec_alpha,
  nChains = 4, nIters = 2000, nCores = 20,
  pathSave = "Results/HARV/Bay_Trait_Q3_correlations.RData",
  Q = 3
)
