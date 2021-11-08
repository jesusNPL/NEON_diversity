###### Part 1 - run taxonomic dimension #####
library(tidyverse)
library(picante)
# remotes::install_github("cran/vegetarian")
library(vegetarian)
library(vegan)

source("R/NEON_diversity/R/Functions/SS_diversity.R")

#### Prepare data #####
load("DATA/matchPhyloComm/matchPhyloComm.RData")
sites <- names(matched_PhyComm)
nSites <- length(sites)

### Spectral species with threshold

commNEON_threshold <- list()
commNEON_observed <- list()

for (i in 1:nSites) {
  print(sites[i])

  # Threshold
  load(paste0(
    "Results/taxoDiversity/SpectralSpecies/SS_",
    sites[i], "_threshold.RData"
  ))

  siteTM_thres <- resThreshold

  siteTM_thres <- siteTM_thres %>%
    mutate(specSpp = paste0("SSpp_", siteTM_thres$SS)) %>%
    select(X, Y, comm, plotID, SS, specSpp, everything())
  # Observed
  load(paste0(
    "Results/taxoDiversity/SpectralSpecies/SS_",
    sites[i], "_observed.RData"
  ))

  siteTM_obs <- resObserved

  siteTM_obs <- siteTM_obs %>%
    mutate(specSpp = paste0("SSpp_", siteTM_obs$SS)) %>%
    select(X, Y, comm, plotID, SS, specSpp, everything())

  ### Make community data matrices for threshold
  commNEON_threshold[[i]] <- make_SSppComm(spectra = siteTM_thres)

  ### Make community data matrices for threshold
  commNEON_observed[[i]] <- make_SSppComm(spectra = siteTM_obs)
}

names(commNEON_threshold) <- sites
names(commNEON_observed) <- sites

### Save spectral communities
save(commNEON_observed,
  file = "Results/taxoDiversity/SpectralCommunities/specComm_NEON_observed.RData"
)

save(commNEON_threshold,
  file = "Results/taxoDiversity/SpectralCommunities/specComm_NEON_threshold.RData"
)

##### Spectral Species diversity under different q #####

#### Prepare data #####
load("DATA/matchPhyloComm/matchPhyloComm.RData")
sites <- names(matched_PhyComm)
nSites <- length(sites)

## Spectral communities observed
load("Results/taxoDiversity/SpectralCommunities/specComm_NEON_observed.RData")
## Spectral communities threshold
load("Results/taxoDiversity/SpectralCommunities/specComm_NEON_threshold.RData")

alphaSpec_thresh_NEON <- list()
alphaSpec_obs_NEON <- list()

for (j in 1:nSites) {

  # Inits
  print(sites[j])

  ## Alpha and beta diversity from community threshold
  commSpecThresh <- commNEON_threshold[[j]]
  plotNames <- rownames(commSpecThresh)
  nPlots <- length(plotNames)

  alphaSpec_thresh_NEON[[j]] <- demon_SSppDIV(
    comm = commSpecThresh,
    q1 = 0, q2 = 1, q3 = 2,
    plotNames = plotNames,
    nPlots = nPlots,
    site = sites[j]
  )

  print("Alpha-beta spectral diversity for threshold, done!!!")

  ## Alpha beta diversity from community observed
  commSpecObs <- commNEON_observed[[j]]
  plotNames <- rownames(commSpecObs)
  nPlots <- length(plotNames)

  alphaSpec_obs_NEON[[j]] <- demon_SSppDIV(
    comm = commSpecObs,
    q1 = 0, q2 = 1, q3 = 2,
    plotNames = plotNames,
    nPlots = nPlots,
    site = sites[j]
  )

  print("Alpha-beta spectral diversity for observed, done!!!")
}

### Aseparate spectral diversity
alphaSpec_thresh_NEON_table <- list()
alphaSpec_obs_NEON_table <- list()

betaSpec_thresh_NEON_table <- list()
betaSpec_obs_NEON_table <- list()

for (i in 1:nSites) {
  alphaSpec_thresh_NEON_table[[i]] <- alphaSpec_thresh_NEON[[i]][[1]]
  alphaSpec_obs_NEON_table[[i]] <- alphaSpec_obs_NEON[[i]][[1]]

  betaSpec_thresh_NEON_table[[i]] <- alphaSpec_thresh_NEON[[i]][[2]]
  betaSpec_obs_NEON_table[[i]] <- alphaSpec_obs_NEON[[i]][[2]]
}

## Combine results
alphaSpec_thresh_NEON_table <- do.call(rbind, alphaSpec_thresh_NEON_table)
alphaSpec_obs_NEON_table <- do.call(rbind, alphaSpec_obs_NEON_table)

betaSpec_thresh_NEON_table <- do.call(rbind, betaSpec_thresh_NEON_table)
betaSpec_obs_NEON_table <- do.call(rbind, betaSpec_obs_NEON_table)

### Save results
save(alphaSpec_obs_NEON_table, alphaSpec_thresh_NEON_table,
  file = "Results/taxoDiversity/alphaSpec_diversity_NEON.RData"
)

save(betaSpec_obs_NEON_table, betaSpec_thresh_NEON_table,
  file = "Results/taxoDiversity/betaSpec_diversity_NEON.RData"
)

##### Alpha and beta diversity from the ground #####

alphaGround_NEON <- list()
betaGround_NEON <- list()

for (j in 1:nSites) {
  print(sites[j])

  ### Inits
  commGround <- matched_PhyComm[[j]][[2]]
  plotNames <- rownames(commGround)
  nPlots <- length(plotNames)

  resSite <- demon_SSppDIV(
    comm = commGround,
    q1 = 0, q2 = 1, q3 = 2,
    plotNames = plotNames,
    nPlots = nPlots,
    site = sites[j]
  )

  alphaGround_NEON[[j]] <- resSite[[1]]
  betaGround_NEON[[j]] <- resSite[[2]]

  print("Alpha-beta diversity for ground, done!!!")
}

alphaGround_NEON_table <- do.call(rbind, alphaGround_NEON)
betaGround_NEON_table <- do.call(rbind, betaGround_NEON)

### Save results
save(alphaGround_NEON_table, alphaGround_NEON,
  file = "Results/taxoDiversity/alphaGround_diversity_NEON.RData"
)

save(betaGround_NEON_table, betaGround_NEON,
  file = "Results/taxoDiversity/betaGround_diversity_NEON.RData"
)
