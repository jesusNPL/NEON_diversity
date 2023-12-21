library(tidyverse)

########### --------------- Taxonomic dimension ---------------- ############

##### Alpha diversity #####

### Spectral species diversity
Spec_diversity <- read_csv("Results/taxoDiversity/Reanalyses/taxo_alphaSS_NEON.csv")

### Ground diversity
Ground_diversity <- read_csv("Results/taxoDiversity/Reanalyses/taxo_alphaGround_NEON.csv")

### Headers
namesSpec <- paste0(names(Spec_diversity[3:ncol(Spec_diversity)]), "_spec")

namesGround <- paste0(names(Ground_diversity[3:ncol(Spec_diversity)]), "_ground")

### Rename columns
names(Spec_diversity) <- c("Site", "plotID", namesSpec) 
head(Spec_diversity)

names(Ground_diversity) <-  c("Site", "plotID", namesGround)
head(Ground_diversity)

### Join diversity metrics from ground and spectra threshold
taxo_spec_thresh_NEON_table <- right_join(
  x = Ground_diversity, 
  y = Spec_diversity, 
  by = c("Site", "plotID")
) 

write_csv(taxo_spec_thresh_NEON_table, 
          file = "Results/RegDATA/Reanalyses/taxo_spec_SS_NEON.csv")

################# ------------------ ##################
##### Combine metrics based spectral distance SAM #####
################ ------------------- ##################

########### --------------- Phylogenetic dimension ---------------- ############

library(tidyverse)

##### Load and arrange data ##### 

## Base data to select plots 
taxo_spec_thresh_NEON_table <- read_csv("Results/RegDATA/Reanalyses/taxo_spec_SS_NEON.csv") 

plot_filter <- paste0(taxo_spec_thresh_NEON_table$Site, 
                      "_", 
                      taxo_spec_thresh_NEON_table$plotID)

### Phylo diversity
load("Results/phyloDiversity/Reanalyses/phylo_DIS_tables.RData")

phylo_NEON_table_q0 

phylo_NEON_table_q1  

phylo_NEON_table_q2 

### Spectral diversity
spec_NEON_table <- read_csv("Results/spectralDiversity/Reanalyses/spectral_normalized_NEON_SAM_q0.csv") 

##### Match tables and prepare regressions #####

##### Spectral diversity
spec_NEON_table <- spec_NEON_table %>%
  mutate(SD2 = SD, MSD2 = MSD) %>%
  dplyr::select(siteID = Site, plotID, SD, SD2, MSD, MSD2, 
                M, mPrime, qHt, qEt, qDT, qDTM)

header_spec <- paste0(names(spec_NEON_table)[3:ncol(spec_NEON_table)], "_spec")

names(spec_NEON_table) <- c("siteID", "plotID", header_spec)

##### Phylogenetic diversity
### Q 0
names(phylo_NEON_table_q0)

## Phylogenetic diversity
phylo_NEON_table_q0 <- phylo_NEON_table_q0 %>%
  dplyr::select(siteID, plotID, PD, PDz, MPD, MPDz, M, mPrime, qHt, qEt, qDT, qDTM)

header_phylo <- paste0(names(phylo_NEON_table_q0)[3:ncol(phylo_NEON_table_q0)], "_phylo")

names(phylo_NEON_table_q0) <- c("siteID", "plotID", header_phylo)

## combine spectral and phylogenetic diversity metrics
phylo_spec_NEON_table_q0 <- full_join(phylo_NEON_table_q0, spec_NEON_table,
                                      by = c("siteID", "plotID")
) %>% 
  mutate(filtro = paste0(siteID, "_", plotID))

### Q 1
phylo_NEON_table_q1 <- phylo_NEON_table_q1 %>%
  dplyr::select(siteID, plotID, PD, PDz, MPD, MPDz, M, mPrime, qHt, qEt, qDT, qDTM)

names(phylo_NEON_table_q1) <- c("siteID", "plotID", header_phylo)

## combine spectral and phylogenetic diversity metrics
phylo_spec_NEON_table_q1 <- full_join(phylo_NEON_table_q1, spec_NEON_table,
                                          by = c("siteID", "plotID")
) %>% 
  mutate(filtro = paste0(siteID, "_", plotID))

### Q 2
phylo_NEON_table_q2 <- phylo_NEON_table_q2 %>%
  dplyr::select(siteID, plotID, PD, PDz, MPD, MPDz, M, mPrime, qHt, qEt, qDT, qDTM)

names(phylo_NEON_table_q2) <- c("siteID", "plotID", header_phylo)

## combine spectral and phylogenetic diversity metrics
phylo_spec_NEON_table_q2 <- full_join(phylo_NEON_table_q2, spec_NEON_table,
                                          by = c("siteID", "plotID")
) %>% 
  mutate(filtro = paste0(siteID, "_", plotID))

### Filter plots based on spectral species for final analyses 
## Q0
phylo_spec_NEON_table_q0 <- phylo_spec_NEON_table_q0 %>% 
  filter(filtro %in% plot_filter) %>% 
  select(!filtro)

## Q1
phylo_spec_NEON_table_q1 <- phylo_spec_NEON_table_q1 %>% 
  filter(filtro %in% plot_filter) %>% 
  select(!filtro)

## Q2
phylo_spec_NEON_table_q2 <- phylo_spec_NEON_table_q2 %>% 
  filter(filtro %in% plot_filter) %>% 
  select(!filtro)

##### Save matched tables ##### 
## Q0
write_csv(phylo_spec_NEON_table_q0, 
          file = "Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q0.csv")

## Q1
write_csv(phylo_spec_NEON_table_q1, 
          file = "Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q1.csv")

## Q2
write_csv(phylo_spec_NEON_table_q2, 
          file = "Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q2.csv")

##### --------------------- Trait diversity METRICS -------------------- ##### 
library(tidyverse)

## Base data to select plots 
taxo_spec_thresh_NEON_table <- read_csv("Results/RegDATA/Reanalyses/taxo_spec_SS_NEON.csv") 

plot_filter <- paste0(taxo_spec_thresh_NEON_table$Site, 
                      "_", 
                      taxo_spec_thresh_NEON_table$plotID) 
##### Match tables and prepare regressions #####

### Spectral diversity
spec_NEON_table <- read_csv("Results/spectralDiversity/Reanalyses/spectral_normalized_NEON_SAM_q0.csv") 

##### Match tables and prepare regressions #####
spec_NEON_table <- spec_NEON_table %>% 
  mutate(SD2 = SD, MSD2 = MSD) %>% 
  select(siteID = Site, plotID, MSD, MSD2, MSm, M, mPrime, qHt, qEt, qDT, qDTM)

header_spec <- paste0(names(spec_NEON_table)[3:ncol(spec_NEON_table)], "_spec")

names(spec_NEON_table) <- c("siteID", "plotID", header_spec)

### Trait diversity
load("Results/traitDiversity/Reanalyses/trait_DIS_tables_Check.RData")

trait_NEON_table_q0 <- trait_NEON_dis_table_q0 %>% 
  select(siteID = Site, plotID, MFD, MFDz, SR, M, mPrime, qHt, qEt, qDT, qDTM)

trait_NEON_table_q1 <- trait_NEON_dis_table_q1 %>% 
  select(siteID = Site, plotID, MFD, MFDz, SR, M, mPrime, qHt, qEt, qDT, qDTM)

trait_NEON_table_q2 <- trait_NEON_dis_table_q2 %>% 
  select(siteID = Site, plotID, MFD, MFDz, SR, M, mPrime, qHt, qEt, qDT, qDTM)

header_trait <- paste0(names(trait_NEON_table_q0)[3:ncol(trait_NEON_table_q0)], "_trait")

names(trait_NEON_table_q0) <- c("siteID", "plotID", header_trait)
names(trait_NEON_table_q1) <- c("siteID", "plotID", header_trait)
names(trait_NEON_table_q2) <- c("siteID", "plotID", header_trait)

### combine spectral and trait diversity metrics

# Q0
trait_spec_NEON_table_q0 <- full_join(x = trait_NEON_table_q0, 
                                      y = spec_NEON_table, 
                                      by = c("siteID", "plotID")) %>% 
  mutate(filtro = paste0(siteID, "_", plotID)) %>% 
  filter(filtro %in% plot_filter) %>% 
  select(!filtro)

# Q1
trait_spec_NEON_table_q1 <- full_join(x = trait_NEON_table_q1, 
                                      y = spec_NEON_table, 
                                      by = c("siteID", "plotID")) %>% 
  mutate(filtro = paste0(siteID, "_", plotID)) %>% 
  filter(filtro %in% plot_filter) %>% 
  select(!filtro)

# Q2
trait_spec_NEON_table_q2 <- full_join(x = trait_NEON_table_q2, 
                                      y = spec_NEON_table, 
                                      by = c("siteID", "plotID")) %>% 
  mutate(filtro = paste0(siteID, "_", plotID)) %>% 
  filter(filtro %in% plot_filter) %>% 
  select(!filtro)

##### Save matched tables ##### 
## Q0
write_csv(trait_spec_NEON_table_q0, 
          file = "Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q0.csv")

## Q1
write_csv(trait_spec_NEON_table_q1, 
          file = "Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q1.csv")

## Q2
write_csv(trait_spec_NEON_table_q2, 
          file = "Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q2.csv")

############### ----------------- ################
##### Combine metrics based spectral species #####
############## ------------------ ################

##### ---------- Spectral dimension based on SS---------- #####
spec_SS_NEON_table <- read_csv("Results/RegDATA/Reanalyses/taxo_spec_SS_NEON.csv") %>% 
  select(siteID = Site, plotID, SRich_spec, Shannon_spec, Simpson_spec)

##### ---------- Phylogenetic dimension ---------- #####

## Q0
phylo_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q0.csv") %>% 
  select(siteID, plotID, qDTM_phylo_q0 = qDTM_phylo)

## Q1
phylo_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q1.csv") %>% 
  select(siteID, plotID, qDTM_phylo_q1 = qDTM_phylo)

## Q2
phylo_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/phylo_spec_DIS_NEON_q2.csv") %>% 
  select(siteID, plotID, qDTM_phylo_q2 = qDTM_phylo)

### Combine phylogenetic diversity metrics
phyloQs <- full_join(x = phylo_NEON_table_q0, 
                     y = phylo_NEON_table_q1, 
                     by = c("siteID", "plotID"))

phyloQs <- full_join(x = phyloQs, 
                     y = phylo_NEON_table_q2, 
                     by = c("siteID", "plotID"))

### Combine phylogenetic and spectral diversity metrics 
phylo_spec_SS_NEON <- full_join(x = phyloQs, 
                             y = spec_SS_NEON_table, 
                             by = c("siteID", "plotID"))

### Save results
write_csv(phylo_spec_SS_NEON, 
          file = "Results/RegDATA/Reanalyses/phylo_spec_SS_NEON.csv")

##### ---------- Trait dimension ---------- #####

## Q0
trait_NEON_table_q0 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q0.csv") %>% 
  select(siteID, plotID, qDTM_trait_q0 = qDTM_trait)

## Q1
trait_NEON_table_q1 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q1.csv") %>% 
  select(siteID, plotID, qDTM_trait_q1 = qDTM_trait)

## Q2
trait_NEON_table_q2 <- read_csv("Results/RegDATA/Reanalyses/trait_spec_DIS_NEON_q2.csv") %>% 
  select(siteID, plotID, qDTM_trait_q2 = qDTM_trait)

### Combine traitgenetic diversity metrics
traitQs <- full_join(x = trait_NEON_table_q0, 
                     y = trait_NEON_table_q1, 
                     by = c("siteID", "plotID"))

traitQs <- full_join(x = traitQs, 
                     y = trait_NEON_table_q2, 
                     by = c("siteID", "plotID"))

### Combine trait and spectral diversity metrics 
trait_spec_SS_NEON <- full_join(x = traitQs, 
                                y = spec_SS_NEON_table, 
                                by = c("siteID", "plotID"))

### Save results
write_csv(trait_spec_SS_NEON, 
          file = "Results/RegDATA/Reanalyses/trait_spec_SS_NEON.csv")
