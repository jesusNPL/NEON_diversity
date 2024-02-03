
############### ------------ GET FIXEF PHYLOGENY  ------------- ############### 

library(tidyverse)

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirPhy <- "Results/Regressions/Reanalyses/phylo_DIS/"

### Inits
## List of RData objects
fits_phy <- list.files(dirPhy)[c(1, 3, 5)]
## Number of RData files
nFiles <- length(fits_phy)
## Diversity orders 
Qs <- c("q0", "q1", "q2")

## List of results
fixef_lst <- list()

for(i in 1:nFiles) {
  
  print(paste0("Getting parameters from ", fits_phy[i], " ...")) 
  
  # Select RData object
  fit_phy <- fits_phy[i]

  # Load results from BLMM
  load(paste0(dirPhy, fit_phy))
  
  # Get fixed effects from BLMM
  fixef <- getFixef(fits = res$fits, 
                    robust = TRUE, 
                    estimates = res$R2_robust,
                    dimension = "phylogeny", 
                    Q = Qs[i], 
                    level = "all") 
  
  # Append results
  fixef_lst[[i]] <- fixef
  
}

### Combine results from 
fixef_table_phy <- do.call(rbind, fixef_lst)

############### ------------ GET FIXEF TRAITS  ------------- ###############

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirTrt <- "Results/Regressions/Reanalyses/trait_DIS/"

### Inits
## List of RData objects
fits_trt <- list.files(dirTrt)[c(1, 3, 5)]
## Number of RData files
nFiles <- length(fits_trt)
## Diversity orders 
Qs <- c("q0", "q1", "q2")

## List of results
fixef_lst_dis <- list()

for(i in 1:nFiles) {
  # Distance based metrics
  print(paste0("Getting parameters from ", fits_trt[i], " ..."))
  
  # Select RData object
  fit_trt_dis <- fits_trt[i] 
  
  # Load results from BLMM
  load(paste0(dirTrt, fit_trt_dis))
  
  # Get fixed effects from BLMM
  fixefdis <- getFixef(fits = res$fits, 
                       robust = TRUE, 
                       estimates = res$R2_robust, 
                       dimension = "trait", 
                       Q = Qs[i], 
                       level = "all") 
  fixef_lst_dis[[i]] <- fixefdis
  
}

fixef_table_trt <- do.call(rbind, fixef_lst_dis)

### Combine distance results
fixef_table_DIS <- bind_rows(fixef_table_phy, fixef_table_trt) %>% 
  mutate(Type = "Distance") %>% 
  select(Type, everything())

############### ------------ GET FIXEF TAXONOMY ------------- ###############

source("R/NEON_diversity/R/Functions/getBEffects.R")

dirTax <- "Results/Regressions/Reanalyses/taxo_SS/"

### Inits
## List of RData objects
fits_tax <- list.files(dirTax)[c(1, 3, 5)]
## Number of RData files
nFiles <- length(fits_tax)
## Diversity orders 
Qs <- c("q0", "q1", "q2")

### Taxonomic dimension SS
load(paste0(dirTax, fits_tax[2]))

fixeftax_SS <- getFixef(fits = res$fits, 
                        robust = TRUE, 
                        estimates = res$R2_robust, 
                        dimension = "taxonomy", 
                        Q = "Q13", 
                        level = "all") 

### Phylogenetic dimension SS
load(paste0(dirTax, fits_tax[1]))

fixefphylo_SS <- getFixef(fits = res$fits, 
                          robust = TRUE, 
                          estimates = res$R2_robust, 
                          dimension = "phylogeny", 
                          Q = "Q13", 
                          level = "all") 

### Trait dimension SS
load(paste0(dirTax, fits_tax[3]))

fixeftrait_SS <- getFixef(fits = res$fits, 
                          robust = TRUE, 
                          estimates = res$R2_robust, 
                          dimension = "trait", 
                          Q = "Q13", 
                          level = "all") 

### Combine results SS
fixef_table_SS <- bind_rows(fixeftax_SS, 
                            fixefphylo_SS, 
                            fixeftrait_SS) %>% 
  mutate(Type = "SS") %>% 
  select(Type, everything())

##### Combine the two data types and save results #####
fixef_ALL <- bind_rows(fixef_table_DIS, 
                       fixef_table_SS)

### Write results 
write_csv(fixef_ALL, 
          file = "Results/Regressions/Reanalyses/FixEffects_table_ALL_results.csv")
