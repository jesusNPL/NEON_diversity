source("R/Functions/SS_diversity.R")

library(dplyr)
library(tidyr)
library(picante)
library(vegetarian)

load("DATA/RData/HARV_SS.RData")

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


ss_th_mean <- ss_threshold[, c(4, 6, 7:432)] %>% 
  group_by(plotID, specSpp) %>% 
  summarise_all("mean")

xx <- ss_threshold %>% 
  count(plotID, specSpp)

head(xx)

##### Make community data matrices ####
harv_comm_thresh <- make_SSppComm(spectra = ss_threshold)

harv_comm_no_thresh <- make_SSppComm(spectra = ss_NO_threshold)

##### SS diversity under different q #####
plotNames <- rownames(harv_comm_thresh)
nPlots <- length(plotNames)

xx <- demon_SSppDIV(comm = harv_comm_thresh, q1 = 0, q2 = 1, q3 = 2, 
                    plotNames = plotNames, nPlots = nPlots)

yy <- demon_SSppDIV(comm = harv_comm_no_thresh, q1 = 0, q2 = 1, q3 = 2, 
                    plotNames = plotNames, nPlots = nPlots)



install.packages("fundiversity")
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

