##### Run functional diversity #####
library(dplyr)
library(tidyr)
library(picante)
library(fundiversity)
library(phytools)
library(gawdis)
source("R/NEON_diversity/R/Functions/AlphaSpecDIV.R")
source("R/NEON_diversity/R/Functions/auxiliar.R")

#### Prepare data #####
load("DATA/RData/HARV_SS.RData") # Spectral species

### Clean vegetation community data 
load("DATA/Harvard_NEON/Data/harv_matchData.RData") # HARV data
rm(harv_taxonomy, matchTrait, matchPhylo)

harv_comm <- harv_comm[harv_comm$plotID_SPEC %in% unique(ss_threshold$plotID), ]
rownames(harv_comm) <- harv_comm$plotID_SPEC

harv_comm_clean <- harv_comm_clean[harv_comm_clean$plotID_SPEC %in% unique(ss_threshold$plotID), ]
rownames(harv_comm_clean) <- harv_comm_clean$plotID_SPEC 

traits <- read.csv("DATA/Traits/Traits_BIEN/GapFilling/imputedTRAITS.csv")
traitHarv <- traits[traits$Species %in% colnames(harv_comm_clean), ]

trtHarv <- data.frame(traitHarv[, c("Species", 
                                    "mean_LAreaPLDryMass", 
                                    "mean_LNCPLDryMass", 
                                    "mean_LDryMass", 
                                    "mean_LCCPLDryMass", 
                                    "mean_LLifeSpan", 
                                    "mean_LCCPLNC", 
                                    "mean_LNCPLArea", 
                                    "mean_LFreshMass")])
# transform to numeric
trtHarv[, 2:9] <- lapply(trtHarv[, 2:9], function(x) as.numeric(as.character(x)))
rownames(trtHarv) <- trtHarv$Species
traitDist <- FD::gowdis(trtHarv[, 2:9])

phy <- ape::read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S3.nex")
phyHarv <- ape::drop.tip(phy, setdiff(phy$tip.label, trtHarv$Species))

harv_comm_match <- match.phylo.comm(phyHarv, harv_comm_clean)

### Run trait diversity 
plotNames <- harv_comm_clean$plotID_SPEC
nPlots <- length(plotNames)

QS <- c(0, 1, 2, 3) # different q values

trait_harv_dis <- list()

for(i in 1:length(QS)) {
  print(paste0("Sarting calculations under Q = ", QS[i]))
  
  trait_harv_dis[[i]] <- demon_TraitDistance(comm = harv_comm_match$comm, 
                                             trait = trtHarv[, 2:9], 
                                             abundance = TRUE, 
                                             Q = QS[i], plotNames = plotNames, 
                                             nPlots = nPlots, gowDist = TRUE, gawDist = FALSE)
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

for(i in 1:length(QS)) {
  
  print(paste0("Starting calculations under ", QS[i], "..."))
  
  spec_harv_dis[[i]] <- demon_DivDistance(spectra = ss_threshold, plotNames = plotNames, 
                                          nPlots = nPlots, specRange = Range, Q = QS[i], 
                                          standardize = FALSE, gowDist = TRUE)
}

spec_harv_table <- do.call(rbind, spec_harv_dis)

### save results 
save(spec_harv_dis, spec_harv_table, file = "Results/HARV/div_spec_harv.RData")

##### Combine trait and spectral distance results #####
Div_trait_spec_alpha <- makeTable_traits(div1 = trait_harv_table, 
                                         div2 = spec_harv_table)

save(Div_trait_spec_alpha, file = "Results/HARV/div_trait_spec_harv.RData")

##### Make Bayesian correlations traits #####
source("R/NEON_diversity/R/Functions/BayesianCorrelation.R")

load("Results/HARV/div_trait_spec_harv.RData")

corr_harv_TRAIT_Q0 <- demon_BayCorrTRAIT(matRES = Div_trait_spec_alpha, 
                                         nChains = 4, nIters = 5000, nCores = 14, 
                                         pathSave = "Results/HARV/Bay_Trait_Q0_correlations.RData", 
                                         Q = 0)

corr_harv_TRAIT_Q1 <- demon_BayCorrTRAIT(matRES = Div_trait_spec_alpha, 
                                         nChains = 4, nIters = 5000, nCores = 14, 
                                         pathSave = "Results/HARV/Bay_Trait_Q1_correlations.RData", 
                                         Q = 1)

corr_harv_TRAIT_Q2 <- demon_BayCorrTRAIT(matRES = Div_trait_spec_alpha, 
                                         nChains = 4, nIters = 5000, nCores = 14, 
                                         pathSave = "Results/HARV/Bay_Trait_Q2_correlations.RData", 
                                         Q = 2)

corr_harv_TRAIT_Q3 <- demon_BayCorrTRAIT(matRES = Div_trait_spec_alpha, 
                                         nChains = 4, nIters = 5000, nCores = 14, 
                                         pathSave = "Results/HARV/Bay_Trait_Q3_correlations.RData", 
                                         Q = 3)
