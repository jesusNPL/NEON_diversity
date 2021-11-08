##### Load required packages #####

packages <- c(
  "fpc", "NbClust", "cluster", "factoextra", "tidyr", "dplyr", "gawdis",
  "svMisc", "parallel", "doParallel"
)

sapply(packages, require, character.only = TRUE)

source("R/NEON_diversity/R/Functions/SpecSpecies.R")

###### Load and prepare data #####
load("DATA/matchPhyloComm/matchPhyloComm.RData")

##### Match spectra wuth community data #####

## Match plots in spectra and ground

specFiles <- list.files(
  path = "DATA/Spectra",
  pattern = "csv"
)

sites <- names(matched_PhyComm)

spec_match <- list()

for (i in 1:length(sites)) {
  print(sites[i])

  specTM <- read.csv(paste0("DATA/Spectra/", sites[i], "_spectra_grid.csv"))

  siteTM <- matched_PhyComm[[i]]$comm

  spec_match[[i]] <- specTM[specTM$plotID %in% rownames(siteTM), ]
}

names(spec_match) <- sites

##### define number of species per plot per site #####
sr_site <- list()

for (i in 1:length(sites)) {
  siteTM <- matched_PhyComm[[i]]$comm

  siteTM[siteTM > 0] <- 1

  srTM <- rowSums(siteTM)
  sr_site[[i]] <- srTM
}

names(sr_site) <- sites

##### Spectral species specification #####

metrics <- c(
  "kl", "ch", "hartigan",
  "cindex", "db", "silhouette",
  # "duda", "pseudot2", "ratkowsky",
  # "beale", "ball",
  "ptbiserial" # , "gap",
  # "frey", "mcclain"#,
  # "gamma", "gplus", "tau", # these are time consuming
  # "hubert",  "dindex", # these are graphical metrics
  #  ccc, scott, marriot, trcovw, tracew, friedman, rubin, # ERROR
  # "dunn", "sdindex", "sdbw"
)

##### Run threshold #####

# resSite <- list()

for (j in 21:length(sites)) {
  print(paste0("Site ", sites[j]))

  ## Inits
  plotTM <- spec_match[[j]]
  plotNames <- unique(plotTM$plotID)
  nPlots <- length(plotNames)

  SR_plots <- sr_site[[j]]

  SR_plots <- SR_plots[names(SR_plots) %in% plotNames]

  SR_plots2 <- SR_plots[order(match(names(SR_plots), plotNames))]

  res <- list()

  for (i in 1:nPlots) {
    print(plotNames[i])
    # for(i in 1:3) {

    ## Filter plot
    specPlot <- plotTM %>%
      filter(plotID == plotNames[i]) # %>%
    # drop_na()

    specPlot <- specPlot[complete.cases(specPlot), ]
    # svMisc::progress(i, max.value = length(plotNames))

    if (nrow(specPlot) <= 100) {
      next
    }

    if (SR_plots2[i] <= 4) {
      SR_plots2[i] <- SR_plots2[i] + 5
    }

    res[[i]] <- demon_ThreshClusters(
      spectra = specPlot[5:430],
      distance = "euclidean",
      Nspp = SR_plots2[i],
      metrics = metrics,
      Ncores = 30
    )

    print(paste0("Thresholds computed for ", plotNames[i], " ... starting new plot"))
  }
  names(res) <- plotNames

  save(res, file = paste0("Results/taxoDiversity/SR_threshold/", sites[j], "_SR_threshold.RData"))
  # resSite[[j]] <- res
}

# SITE GRSM THROUGHOUT ERRORS

##### extract thresholds #####

thrs_NEON <- list()

for (j in 1:length(sites)) {
  print(paste0("Site ", sites[j]))

  load(paste0("Results/taxoDiversity/SR_threshold/", sites[j], "_SR_threshold.RData"))

  plotNames <- names(res)

  nPlots <- length(plotNames)

  thrs <- numeric(length = nPlots)

  for (i in 1:nPlots) {
    pt_th <- res[[i]]

    if (is.null(pt_th) == TRUE) {
      next
    }

    thrs[i] <- pt_th$Mean_rule
  }

  names(thrs) <- plotNames

  thrs_NEON[[j]] <- thrs

  Sys.sleep(3)
}

names(thrs_NEON) <- sites

### Match SR with the plots in the observed community ###

sr_site_NEON <- list()

for (k in 1:length(sites)) {
  print(paste0("Site ", sites[k]))

  srTHR <- thrs_NEON[[k]]

  srTM <- sr_site[[k]]

  SR_plots <- srTM[names(srTM) %in% names(srTHR)]

  SR_plots2 <- SR_plots[order(match(names(SR_plots), names(srTHR)))]

  sr_site_NEON[[k]] <- SR_plots2
}

names(sr_site_NEON) <- sites

save(thrs_NEON, sr_site_NEON,
  file = "Results/taxoDiversity/SR_plots_NEON.RData"
)

##### Run SPECTRAL SPECIES #####
library(tidyverse)
library(cluster)

source("R/NEON_diversity/R/Functions/SpecSpecies.R")

### Load and prepare data
load("DATA/matchPhyloComm/matchPhyloComm.RData")
sites <- names(matched_PhyComm)

load("Results/taxoDiversity/SR_plots_NEON.RData")

ss_threshold <- list()
ss_observed <- list()

for (i in 1:length(sites)) {
  print(paste0("Site ", sites[i]))

  srOBS <- sr_site_NEON[[i]]
  srTHR <- thrs_NEON[[i]]

  plotNames <- names(srTHR)

  nPlots <- length(plotNames)

  specTM <- read.csv(paste0("DATA/Spectra/", sites[i], "_spectra_grid.csv"))

  specTM_clean <- specTM %>%
    filter(plotID %in% plotNames)

  ## Run SS with thresholds
  resThreshold <- demon_SS(
    spectra = specTM_clean, # spectra data
    threshold = srTHR, # threshold, i.e., number of clusters per plot
    nPlots = nPlots, # number of plots
    plotNames = names(srTHR)
  ) # plot names

  save(resThreshold, file = paste0(
    "Results/taxoDiversity/SpectralSpecies/SS_",
    sites[i], "_threshold.RData"
  ))

  ss_threshold[[i]] <- resThreshold
  ## Run SS with the number of species in the observed communities
  resObserved <- demon_SS(
    spectra = specTM_clean,
    threshold = srOBS, # number of species per plot
    nPlots = nPlots,
    plotNames = names(srOBS)
  )

  save(resObserved, file = paste0(
    "Results/taxoDiversity/SpectralSpecies/SS_",
    sites[i], "_observed.RData"
  ))

  ss_observed[[i]] <- resObserved
}

## Save spectral species
names(ss_threshold) <- sites
names(ss_observed) <- sites

save(ss_threshold, ss_observed,
  file = "Results/taxoDiversity/SpectralSpecies/ALL_sites_SS_NEON.RData"
)
