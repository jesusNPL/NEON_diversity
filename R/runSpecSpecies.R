##### Load required packages #####

packages <- c("fpc", "NbClust", "cluster", "factoextra", "tidyr", "dplyr", "gawdis", 
              "svMisc", "parallel", "doParallel")

sapply(packages, require, character.only = TRUE)

source("R/NEON_diversity/R/Functions/SpecSpecies.R")

###### Load and prepare data #####

harv_veg <- read.csv("DATA/Harvard_NEON/Data/HARV_veg_cover_compiled.csv")

harv_spec <- read.csv("DATA/Harvard_NEON/Data/HARV_spectra_grid.csv")

## Match plots in spectra and ground
harv_spec <- harv_spec[harv_spec$plotID %in% harv_veg$plotID_SPEC, ]

names(harv_spec)

# define number of species per plot
PAM <- harv_veg[, 6:146]
PAM[PAM > 0 ] <- 1

SR_plots <- rowSums(PAM)
names(SR_plots) <- harv_veg$plotID_SPEC

##### Spectral species specification #####

metrics <- c("kl", "ch", "hartigan", 
             "cindex", "db", "silhouette", 
             #"duda", "pseudot2", "ratkowsky", 
             #"beale", "ball", 
             "ptbiserial"#, "gap", 
             #"frey", "mcclain"#, 
             # "gamma", "gplus", "tau", # these are time consuming 
             #"hubert",  "dindex", # these are graphical metrics
             #  ccc, scott, marriot, trcovw, tracew, friedman, rubin, # ERROR
             #"dunn", "sdindex", "sdbw"
             )

plotNames <- harv_veg$plotID_SPEC

nPlots <- length(plotNames)

##### Run threshold #####
res <- list()

for(i in 1:nPlots) {
  #for(i in 1:3) {
    
  print(plotNames[i])
  #svMisc::progress(i, max.value = length(plotNames)) 
  
  if (plotNames[i] == "P_23") {
    next
  }
  
  specPlot <- harv_spec %>% filter(plotID == plotNames[i])
  

  res[[i]] <- demon_ThreshClusters(spectra = na.omit(specPlot[5:430]), 
                                   distance = "euclidean", 
                                   Nspp = SR_plots[i], 
                                   metrics = metrics, 
                                   Ncores = 30)
  
  print(paste0("Thresholds computed for ", plotNames[i], " ... starting new plot"))
  
} 

names(res) <- plotNames

## remove NULL objects
Filter(Negate(is.null), res)

### Save thresholds #####

save(res, file = "DATA/RData/HARV_thresholds.RData")

##### Run spectral species #####

thrs <- numeric(length = nPlots)

for(i in 1:nPlots) { 
  
  pt_th <- res[[i]]
  
  if (is.null(pt_th) == TRUE) {
    next
  }
  
  thrs[i] <- pt_th$Mean_rule
}

names(thrs) <- plotNames

thrs
## remove plotID = P_23 because the large amount of NA

harv_spec_clean <- harv_spec %>% 
  filter(plotID != "P_23")

thrs_2 <- thrs[-17]
thrs_2

SR_plots_2 <- SR_plots[-17]
SR_plots_2

## Run SS with thresholds
ss_threshold <- demon_SS(spectra = harv_spec_clean, # spectra data
                         threshold = thrs_2, # threshold, i.e., number of clusters per plot
                         nPlots = nPlots-1, # number of plots
                         plotNames = names(thrs_2)) # plot names

## Run SS with the number of species in the observed communities
ss_NO_threshold <- demon_SS(spectra = harv_spec_clean, 
                            threshold = SR_plots_2, # number of species per plot
                            nPlots = nPlots-1, 
                            plotNames = names(SR_plots_2))

## Save spectral species
save(ss_NO_threshold, ss_threshold, 
     file = "Data/RData/HARV_SS.RData")

