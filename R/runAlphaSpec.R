
library(dplyr)
library(tidyr)
library(picante)
#remotes::install_github("cran/vegetarian")
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

spec_HARV_div_thresh <- demon_SSppDIV(comm = harv_commSPEC_thresh, q1 = 0, q2 = 1, q3 = 2, 
                    plotNames = plotNames, nPlots = nPlots)

spec_HARV_div_no_thresh <- demon_SSppDIV(comm = harv_commSPEC_no_thresh, q1 = 0, q2 = 1, q3 = 2, 
                    plotNames = plotNames, nPlots = nPlots)

##### Vegetation diversity under different q ##### 

veg_HARV_div <- demon_SSppDIV(comm = harv_comm_clean[, 4:121], q1 = 0, q2 = 1, q3 = 2, 
                                      plotNames = plotNames, nPlots = nPlots)

##### Combine information  of the two dimensions #####
Div_tax_alpha_thresh <- makeTable_taxonomy(div1 = veg_HARV_div$commDiv, 
                                   div2 = spec_HARV_div_thresh$commDiv)

Div_tax_alpha_thresh$divOBS$herv_cov <- harv_comm$herb_cov
Div_tax_alpha_thresh$divOBS$tree_cov <- harv_comm$tree_cov

Div_tax_alpha_no_thresh <- makeTable_taxonomy(div1 = veg_HARV_div$commDiv, 
                                      div2 = spec_HARV_div_no_thresh$commDiv)

Div_tax_alpha_no_thresh$divOBS$herv_cov <- harv_comm$herb_cov
Div_tax_alpha_no_thresh$divOBS$tree_cov <- harv_comm$tree_cov

### Save tables with alpha calculations
saveRDS(Div_tax_alpha_no_thresh, file = "Results/HARV/taxo_alpha_no_thresh.rds")
saveRDS(Div_tax_alpha_thresh, file = "Results/HARV/taxo_alpha_thresh.rds")

save(Div_tax_alpha_thresh, Div_tax_alpha_no_thresh, 
     file = "Results/HARV/div_tax_harv.RData")

##### Make Bayesian robust correlations #####
source("R/NEON_diversity/R/Functions/BayesianCorrelation.R")

load("Results/HARV/div_tax_harv.RData")

corr_harv_tax_thresh <- demon_BayCorrTAX(matRES = Div_tax_alpha_thresh$divOBS, 
                                         nChains = 4, nIters = 2000, nCores = 20, 
                                         pathSave = "Results/HARV/Bay_TAX_thresh_correlations.RData")

corr_harv_tax_no_thresh <- demon_BayCorrTAX(matRES = Div_tax_alpha_no_thresh$divOBS, 
                                            nChains = 4, nIters = 2000, nCores = 20, 
                                            pathSave = "Results/HARV/Bay_TAX_no_thresh_correlations.RData")






library(vegan)
data(BCI)
i <- sample(nrow(BCI), 12)
mod <- renyi(BCI[i,])
plot(mod)
mod <- renyiaccum(BCI[i,])
plot(mod, as.table=TRUE, col = c(1, 2, 2))
persp(mod)


data(BCI)
H <- diversity(BCI)
simp <- diversity(BCI, "simpson")
invsimp <- diversity(BCI, "inv")
## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
unbias.simp <- rarefy(BCI, 2) - 1
## Fisher alpha
alpha <- fisher.alpha(BCI)
## Plot all
pairs(cbind(H, simp, invsimp, unbias.simp, alpha), pch="+", col="blue")
## Species richness (S) and Pielou's evenness (J):
S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
J <- H/log(S)
## beta diversity defined as gamma/alpha - 1:
data(dune)
data(dune.env)
alpha <- with(dune.env, tapply(specnumber(dune), Management, mean))
gamma <- with(dune.env, specnumber(dune, Management))
gamma/alpha - 1


#install.packages("fundiversity")
library(fundiversity)

data(traits_birds)
fd_raoq(traits_birds)
fd_fric(traits_birds)
fd_feve(traits_birds)
fd_fdiv(traits_birds)

data("traits_plants")
fd_raoq(traits_plants)
fd_fric(traits_plants)
fd_feve(traits_plants)
fd_fdiv(traits_plants)


remotes::install_github("ibartomeus/fundiv")
library(fundiv)


ex1 <- betaFD(c1 = c("sp3", "sp2", "sp1", "sp4", "sp5"), 
              c2 = c("sp6", "sp7", "sp8", "sp4", "sp5"), 
              S = dummy$trait)
ex1

fd <- FD_dendro2(S, A, Tree, Cluster.method = "average", 
                 ord = "podani")

ex1 <- FD_dendro2(A = dummy$abun, Tree = FDtree(S = dummy$trait, w = NA,
                                                Distance.method = "gower", 
                                                ord = "podani", 
                                                Cluster.method = "average"))
ex1

Xtree(h)
ex1 <- FD_Clark(S = dummy$trait, A = dummy$abun, Cluster.method = "average", ord = "podani",
                Weigthedby = "abundance")
ex1


Tree = FDtree(S = dummy$trait, w = NA,
              Distance.method = "gower", 
              ord = "podani", 
              Cluster.method = "average")

Xtree(Tree$h2.prime)


Ea0(A = dummy$abun)
Ea0(A = dummy$abun, scales = c(0.25,0.5,1,2,4,8,Inf))
EinvD(A = dummy$abun)
EJ(A = dummy$abun)
Evar(A = dummy$abun)
#calculate all of them
eve1 <- Eve(A = dummy$abun)



ss_th_mean <- ss_threshold[, c(4, 6, 7:432)] %>% 
  group_by(plotID, specSpp) %>% 
  summarise_all("mean")

ss_pt1 <- ss_th_mean %>% filter(plotID == "P_1")

ss <- data.frame(ss_pt1[, 3:428])
rownames(ss) <- ss_pt1$specSpp

dis <- gowdis(ss)
tree <- hclust(dis, method = "average")
tree <- as.phylo(tree)
plot(tree)


setwd("../Dropbox/Macrosystems_NEON/")

load("DATA/RData/HARV_SS.RData")
library(tidyr)
library(dplyr)

harv_p1 <- ss_NO_threshold %>% 
  filter(plotID == "P_1")

harv_p1 <- harv_p1[, c(1, 2, 5:430)]
harv_p1[1:10, 1:7]

xx <- demon_SpecBeta(spectra = harv_p1, coords = harv_p1[, 1:2], 
                     window = 1.5, metric = "bray")

raoPICANTE <- raoD(harv_p1[, 5:430])
eveFUNDIV <- Eve(harv_p1[, 5:430])
Sys.time()

xx$specBeta
xx$specBetaPoints
xx$specBetaRaster

plot(xx$specBetaRaster$specTurnover)
plot(xx$specBetaRaster$specNestedness)
plot(xx$specBetaRaster$specBetaRatio)
plot(xx$specBetaRaster$specBeta)

plot(xx$specBetaRaster)


xx <- demon_specFour(spectra = harv_p1, coords = harv_p1[, 1:2], window = 2, nPCA = 15)


pts <- unique(ss_NO_threshold$plotID)
npts <- length(pts)


xx
varcompHeight <- ape::varcomp(xx, scale = 1)
varcompHeight

setwd("C:/Users/jpintole/Dropbox/Macrosystems_NEON")
load("DATA/RData/HARV_SS.RData")
library(tidyr)
library(dplyr)
require(vegan)

#devtools::install_github("ibartomeus/fundiv")

source("R/Functions/AlphaSpecDIV.R")

xx <- demon_specFour(spectra = harv_p1, coords = harv_p1[, 1:2], window = 2, nPCA = 15)


pts <- unique(ss_NO_threshold$plotID)
npts <- length(pts)

xxx <- demon_DivDistance(spectra = ss_NO_threshold, distance = NULL, 
                         standardize = TRUE, gowDist = TRUE, Q = 1, 
                         nPlots = npts, 
                         plotNames = pts, specRange = 5:430) 
xxx

yyy <- demon_SpecRaoEve(spectra = ss_NO_threshold, standardize = FALSE, 
                             nPlots = npts, plotNames = pts, specRange = 5:430)

yyy

RAO_EVE_HARV <- list()

for(i in 1:npts) {
  
  print(plotNames[i])
  
  harv <- ss_NO_threshold %>% filter(plotID == plotNames[i])
  
  RAO_EVE_HARV[[i]] <- demon_SpecRaoEve_RASTER(spectra = harv, coords = 1:2, 
                                               specRange = 5:430, standardize = FALSE, 
                                               plotName = plotNames[i])
  
  print(paste0("rasters calculated for ", plotNames[i], " starting new plot..."))
}

