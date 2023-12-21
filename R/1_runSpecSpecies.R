##### Load required packages #####

packages <- c(
  "fpc", "NbClust", "cluster", 
  "factoextra", "tidyr", "dplyr",
  "svMisc", "parallel", "doParallel"
)

sapply(packages, require, character.only = TRUE)

source("R/NEON_diversity/R/Functions/SpecSpecies.R")

###### ---------- Step 1 - Load and prepare data ---------- #####
load("DATA/matchPhyloComm/matchPhyloComm.RData")

##### Match spectra with community data #####

## Match plots in spectra and ground

specFiles <- list.files(
  path = "DATA/Spectra/normalized",
  pattern = "csv"
)

sites <- names(matched_PhyComm)

spec_match <- list()

for (i in 1:length(sites)) { 
  
  print(sites[i])

  # spectra data normalized
  specTM <- read.csv(paste0("DATA/Spectra/normalized/normalized_mask_", sites[i], "_spectra_grid.csv"))

  siteTM <- matched_PhyComm[[i]]$comm

  spec_match[[i]] <- specTM[specTM$plotID %in% rownames(siteTM), ] 
  
}

names(spec_match) <- sites

##### define number of species per plot per site ##### 
# Information to be used in setting maximum number of potential spectral clusters
sr_site <- list()

for (i in 1:length(sites)) { 
  
  siteTM <- matched_PhyComm[[i]]$comm

  siteTM[siteTM > 0] <- 1

  srTM <- rowSums(siteTM)
  
  sr_site[[i]] <- srTM 
  
}

names(sr_site) <- sites

save(sr_site, 
     file = "DATA/Spectral_distance_SAM/SR_obs_plot_by_Site.RData")

##### ---------- Step 2 - Spectral species specification ---------- #####

### Select metrics to estimate the potential number of clusters per plot 
metrics <- c(
  "kl", "ch", "hartigan",
  "cindex", "db", "silhouette",
  "duda", "pseudot2", "ratkowsky",
  "beale", "ball",
  "ptbiserial", #"gap",
  "frey", "mcclain",
  "gamma", "gplus", "tau", # these are time consuming
  # "hubert",  "dindex", # these are graphical metrics
  #ccc, scott, marriot, trcovw, tracew, friedman, rubin, # ERROR
  "dunn", "sdindex", "sdbw"
)

### Selected clustering metrics
metricsTen <- c(
  "kl", "ch", "hartigan", "cindex", "db", 
  "silhouette", "ptbiserial", "duda", "tau", "dun"
)

##### Run threshold ##### 
library(tidyverse) # Easily Install and Load the 'Tidyverse', CRAN v2.0.0 

# Load number of species per plot by site
load("DATA/Spectral_distance_SAM/SR_obs_plot_by_Site.RData")

### Inits

# Site names
siteNames <- sites
# Number of sites
nSites <- length(siteNames)
# number of metrics
nMetrics <- length(metricsTen)

### Start calculations
for (j in 1:nSites) { 
  
  print(paste0("Site ", siteNames[j]))

  ## Inits
  plotTM <- spec_match[[j]]
  plotNames <- unique(plotTM$plotID)
  nPlots <- length(plotNames) 
  
  ## Match plot names in both spectra and richness
  SR_plots <- sr_site[[j]]

  SR_plots <- SR_plots[names(SR_plots) %in% plotNames] 
  
  # Match names and positions between plots 
  SR_plots3 <- SR_plots[order(match(names(SR_plots), plotNames))] # keep original observed SPP for sanity check!
  SR_plots2 <- SR_plots[order(match(names(SR_plots), plotNames))] # after using median SR for plots with less than 4 species
  
  # Store temporary results
  res <- list() 
  
  # Start computations within plots by site
  for (i in 1:nPlots) { 
    
    print(plotNames[i])
    
    ## Filter plot
    specPlot <- plotTM %>% 
      filter(plotID == plotNames[i])
    
    ## remove NAs that correspond to masked pixels  
    specPlot <- specPlot %>% 
      drop_na() 
    
    ## if the spectra data contain less than 10 rows (or 10 x 10 meters) skip 
    if (nrow(specPlot) <= 10) { 
      
      next 
      
    } 
    
    ## If the plot has less than 4 species add the MEDIAN of richness across the site to allow the computation of maximum potential clusters 
    # for example, the median SR for site OSBS is 15 (mean = 21) and the standard deviation is 13
    if (SR_plots2[i] <= 4) { 
      
      SR_plots2[i] <- round(median(SR_plots2), 0) # adding median SR across the site to plots with less than 4 species
      
    } 
    
    # sample.int(m, k) : cannot take a sample larger than the population when 'replace = FALSE'
    if (nrow(specPlot) < SR_plots2[i]) { 
      
      next 
      
    }

    ## Run aux function that allow the computation of multiple cluster metrics 
    resTM <- demon_ThreshClusters( 
      spectra = specPlot[4:ncol(specPlot)], # spectra data after NDVI and NIR masks
      distance = "euclidean", # distance 
      Nspp = SR_plots2[i], # maximum number of species based on the observed richness in the plot
      metrics = metricsTen, # vector with clustering metric names
      nMetrics = nMetrics, # number of clustering metrics
      Ncores = 20
    ) 
    
    # Add supplementary information to the resulting data.frame
    resTM$N_spp_orig <- SR_plots3[i] # adding original species richness 
    resTM$N_spp_with_median <- SR_plots2[i] # adding median SR "IF" the plot has less than 4 species
    resTM$plotID <- plotNames[i] 
    
    res[[i]] <- resTM

    print(paste0("Thresholds computed for ", plotNames[i], " ... starting new plot")) 
    
  } # end for i
  
  results <- do.call(rbind, res)
  
  results$siteID <- siteNames[j] # Assign site names
  
  # write results by site
  write_csv(results, 
            file = paste0("Results/taxoDiversity/Reanalyses/thresholds/", sites[j], "_threshold_ten.csv"))

} # end for j

##### ---------- Step 3 - Threshold calculation or expected number of clusters ---------- #####

### Calculate mean number of clusters
make_rules <- function(resCluster, site) { 
  
  ## Descriptive rules
  res1 <- resCluster %>% 
    group_by(plotID) %>% 
    summarize(mean_rule = round(mean(N_clusters), 0), 
              median_rule = round(median(N_clusters), 0), 
              sd_rule = round(sd(N_clusters), 0), 
              max_rule = max(N_clusters), 
              min_rule = min(N_clusters)) 
  
  ## Majority consensus rule
  res2 <- resCluster %>% 
    group_by(plotID) %>% 
    count(N_clusters) %>% 
    summarise(majority_rule = max(n)) 
  
  ## Combine results
  thresholds <- full_join(res1, res2, by = "plotID") %>% 
    mutate(siteID = site) %>% 
    select(siteID, everything())
  
  return(thresholds)
  
}

### Inits
# Number of clusters per metric - using ten metrics
fileNames <- list.files(
  path = "Results/taxoDiversity/Reanalyses/thresholds/", 
  pattern = "csv" 
  )

# Site names 
siteNames <- fileNames %>% 
  tibble(file = .) %>% 
  separate(file, c("siteID"))

siteNames <- siteNames$siteID

# Number of sites
nSites <- length(siteNames)

## Global rules
rules <- list()

### Run rules 
for(i in 1:nSites) { 
  
  print(siteNames[i])
  
  ## Load thresholds by metric
  dt <- read_csv(paste0("Results/taxoDiversity/Reanalyses/thresholds/", 
                        fileNames[i])) 
  
  ## Select seven metrics as the original analyses
  dt <- dt %>% 
    filter(Metric %in% c("kl", "ch", "hartigan", "cindex", 
                         "db", "silhouette", "ptbiserial")) 
  
  ## Run threshold rules
  res <- make_rules(resCluster = dt, 
                    site = siteNames[i])
  
  rules[[i]] <- res 
  
  ## Save results per site
  write_csv(res, 
            file = paste0("Results/taxoDiversity/SR_threshold/Reanalyses/", 
                          siteNames[i], "_7metrics_rules.csv"))
  
}

### Combine results 
rules_NEON <- do.call(rbind, rules)

### Save results
write_csv(rules_NEON, 
          file = "Results/taxoDiversity/SR_threshold/Reanalyses/NEON_7metrics_rules.csv")

### Explore results using both threshold specifications - 10 vs 7 metrics
m7 <- read_csv("Results/taxoDiversity/SR_threshold/Reanalyses/NEON_7metrics_rules.csv")
m10 <- read_csv("Results/taxoDiversity/SR_threshold/Reanalyses/NEON_10metrics_rules.csv")

cor.test(m7$mean_rule, m10$mean_rule)
plot(m7$mean_rule, m10$mean_rule)

cor.test(m7$max_rule, m10$max_rule)
plot(m7$max_rule, m10$max_rule)

##### ---------- Step 4 - Create spectral species using thresholds ---------- #####

library(tidyverse) # Easily Install and Load the 'Tidyverse', CRAN v2.0.0 
library(cluster) 

### Auxiliary functions
source("R/NEON_diversity/R/Functions/SpecSpecies.R")

### Load and prepare data
load("DATA/matchPhyloComm/matchPhyloComm.RData")

## Site names
sites <- names(matched_PhyComm)

## Thresholds by site and plots
threshold10 <- read_csv("Results/taxoDiversity/SR_threshold/Reanalyses/NEON_10metrics_rules.csv")

### Run spectral species
for (i in 1:length(sites)) { 
  
  print(sites[i])
  
  ### Threshold by site
  thresh_site <- threshold10 %>% 
    filter(siteID == sites[i])
  
  ### spectra data normalized
  specTM <- read.csv(paste0("DATA/Spectra/normalized/normalized_mask_", sites[i], "_spectra_grid.csv"))
  
  # Filter plots that contain thresholds
  specTM <- specTM %>% 
    filter(plotID %in% thresh_site$plotID)
  
  ### Inits
  # Plot names
  plotNames <- unique(specTM$plotID)
  # Number of plots
  nPlots <- length(plotNames)
  
  ### Match thresholds ID with plotID 
  thresh_plot <- setNames(thresh_site$mean_rule, thresh_site$plotID)
  
  thresh_plot <- thresh_plot[order(match(names(thresh_plot), plotNames))]
  
  ### Run SS with thresholds
  resThreshold <- demon_SS(
    spectra = specTM, # spectra data
    threshold = thresh_plot, # threshold, i.e., number of clusters per plot
    nPlots = nPlots, # number of plots
    plotNames = names(thresh_plot) # plot names
  ) 
  
  write_csv(resThreshold, 
            file = paste0(
              "Results/taxoDiversity/SpectralSpecies/Reanalyses/SS_", 
              sites[i], "_threshold_spectra.csv"
  ) )
  
}

##### ---------- Step 5 - Create spectral communities ---------- #####

library(tidyverse) # Easily Install and Load the 'Tidyverse', CRAN v2.0.0 

### Auxiliary functions
source("R/NEON_diversity/R/Functions/SS_diversity.R")

### Load and prepare data
load("DATA/matchPhyloComm/matchPhyloComm.RData")

## Site names
siteNames <- names(matched_PhyComm)

## Number of sites
nSites <- length(siteNames)

### Run spectral communities
for(j in 1:nSites) { 
  
  print(siteNames[j]) 
  
  ### Load spectra with spectral species
  spec_SS <- read_csv(paste0("Results/taxoDiversity/SpectralSpecies/Reanalyses/SS_", 
                        siteNames[j], "_threshold_spectra.csv")) 
  
  ### Create community data matrix using spectral species
  comm_SS <- make_SSppComm(SS_spectra = spec_SS)
  
  ### Save communities by site
  write_csv(
    x = comm_SS, 
    file = paste0("Results/taxoDiversity/SpectralCommunities/Reanalyses/spectral_CDM_", 
                  siteNames[j], "_10metrics.csv"))
  
  print(paste0("CDM based on spectral species for site ", siteNames[j], " created successfully!"))
  
}

