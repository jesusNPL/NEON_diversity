library(tidyverse)

##### Load and arrange data #####
### Phylo diversity
load("Results/phyloDiversity/phylo_DIS_tables.RData")

phylo_NEON_table_q0 <- phylo_NEON_table_q0 %>%
  select(Site, plotID, everything())

phylo_NEON_table_q1 <- phylo_NEON_table_q1 %>%
  select(Site, plotID, everything())

phylo_NEON_table_q2 <- phylo_NEON_table_q2 %>%
  select(Site, plotID, everything())

phylo_NEON_table_q3 <- phylo_NEON_table_q3 %>%
  select(Site, plotID, everything())

### Spectral diversity
spec_NEON_div <- readRDS("Results/spectralDiversity/spectral_NEON_q0.rds")

spec_NEON_div <- spec_NEON_div %>%
  select(Site, plotID, everything())

##### Match tables and prepare regressions #####

##### Spectral diversity
spec_NEON_div_sel <- spec_NEON_div %>%
  select(Site, plotID, SD, MSD, M, mPrime, qHt, qEt, qDT, qDTM, MSm, Range) %>%
  mutate(SD2 = SD, MSD2 = MSD) %>%
  select(Site, plotID, SD, SD2, MSD, MSD2, everything())

names(spec_NEON_div_sel) <- c(
  "Site", "plotID", "SD", "SD2", "MSD", "MSD2",
  "M_s", "mPrime_s", "qHt_s", "qEt_s", "qDT_s", "qDTM_s",
  "MSm", "Range"
)

##### Phylogenetic diversity
### Q 0
names(phylo_NEON_table_q0)
names(spec_NEON_div)

## Phylogenetic diversity
phylo_NEON_table_q0_sel <- phylo_NEON_table_q0 %>%
  select(Site, plotID, PD, PDz, MPD, MPDz, M, mPrime, qHt, qEt, qDT, qDTM)

names(phylo_NEON_table_q0_sel) <- c(
  "Site", "plotID", "PD", "PDz", "MPD", "MPDz",
  "M_p", "mPrime_p", "qHt_p", "qEt_p", "qDT_p", "qDTM_p"
)

## combine spectral and phylogenetic diversity metrics
phylo_spec_NEON_table_q0 <- full_join(phylo_NEON_table_q0_sel, spec_NEON_div_sel,
  by = c("Site", "plotID")
)

### Q 1
names(phylo_NEON_table_q1)
names(spec_NEON_div)

phylo_NEON_table_q1_sel <- phylo_NEON_table_q1 %>%
  select(Site, plotID, PD, PDz, MPD, MPDz, M, mPrime, qHt, qEt, qDT, qDTM)

names(phylo_NEON_table_q1_sel) <- c(
  "Site", "plotID", "PD", "PDz", "MPD", "MPDz",
  "M_p", "mPrime_p", "qHt_p", "qEt_p", "qDT_p", "qDTM_p"
)

## combine spectral and phylogenetic diversity metrics
phylo_spec_NEON_table_q1 <- full_join(phylo_NEON_table_q1_sel, spec_NEON_div_sel,
  by = c("Site", "plotID")
)

### Q 2
names(phylo_NEON_table_q2)
names(spec_NEON_div)

phylo_NEON_table_q2_sel <- phylo_NEON_table_q2 %>%
  select(Site, plotID, PD, PDz, MPD, MPDz, M, mPrime, qHt, qEt, qDT, qDTM)

names(phylo_NEON_table_q2_sel) <- c(
  "Site", "plotID", "PD", "PDz", "MPD", "MPDz",
  "M_p", "mPrime_p", "qHt_p", "qEt_p", "qDT_p", "qDTM_p"
)

## combine spectral and phylogenetic diversity metrics
phylo_spec_NEON_table_q2 <- full_join(phylo_NEON_table_q2_sel, spec_NEON_div_sel,
  by = c("Site", "plotID")
)

### Q 3
names(phylo_NEON_table_q3)
names(spec_NEON_div)

phylo_NEON_table_q3_sel <- phylo_NEON_table_q3 %>%
  select(Site, plotID, PD, PDz, MPD, MPDz, M, mPrime, qHt, qEt, qDT, qDTM)

names(phylo_NEON_table_q3_sel) <- c(
  "Site", "plotID", "PD", "PDz", "MPD", "MPDz",
  "M_p", "mPrime_p", "qHt_p", "qEt_p", "qDT_p", "qDTM_p"
)

## combine spectral and phylogenetic diversity metrics
phylo_spec_NEON_table_q3 <- full_join(phylo_NEON_table_q3_sel, spec_NEON_div_sel,
  by = c("Site", "plotID")
)

##### Save matched tables #####

save(phylo_spec_NEON_table_q0, phylo_spec_NEON_table_q1,
  phylo_spec_NEON_table_q2, phylo_spec_NEON_table_q3,
  file = "Results/RegDATA/phylo_spec_DIS_NEON.RData"
)

##### ------------------------- MOMENTS ---------------------- #####
library(tidyverse)

##### Match tables moments ######

### Load data
load("Results/moments/moments_NEON.RData")

moments_Spec_neon <- readRDS("Results/spectralDiversity/spectral_NEON_moments.rds")

### Combine data
nSites <- length(moments_Spec_neon)

momentTaxoSpec <- list()
momentTraitSpec <- list()
momentPhyloSpec <- list()

for (i in 1:nSites) {
  # Taxonomy
  momentTaxoSpec[[i]] <- full_join(moments_NEON[[i]][[1]], moments_Spec_neon[[i]],
    by = c("Site", "plotID")
  )
  # Traits
  momentTraitSpec[[i]] <- full_join(moments_NEON[[i]][[2]], moments_Spec_neon[[i]],
    by = c("Site", "plotID")
  )
  # Phylogeny
  momentPhyloSpec[[i]] <- full_join(moments_NEON[[i]][[3]], moments_Spec_neon[[i]],
    by = c("Site", "plotID")
  )
}

momentTaxoSpec_table <- do.call(rbind, momentTaxoSpec)
momentTraitSpec_table <- do.call(rbind, momentTraitSpec)
momentPhyloSpec_table <- do.call(rbind, momentPhyloSpec)

### Combine tables moments
tbl_moments_NEON <- list(
  taxoSpec = momentTaxoSpec_table,
  traitSpec = momentTraitSpec_table,
  phyloSpec = momentPhyloSpec_table
)

##### Save matched tables #####

save(tbl_moments_NEON,
  file = "Results/RegDATA/moments_ALL_NEON.RData"
)

##### ------------------- Trait diversity ----------------------- #####

library(tidyverse)

### Spectral diversity
spec_NEON_div <- readRDS("Results/spectralDiversity/spectral_NEON_q0.rds")

spec_NEON_div <- spec_NEON_div %>%
  select(Site, plotID, everything())

##### Match tables and prepare regressions #####

##### Spectral diversity
spec_NEON_div_sel <- spec_NEON_div %>%
  select(Site, plotID, SD, MSD, M, mPrime, qHt, qEt, qDT, qDTM, MSm, Range) %>%
  mutate(SD2 = SD, MSD2 = MSD) %>%
  select(Site, plotID, MSD, MSD2, MSm, M, mPrime, qHt, qEt, qDT, qDTM)

names(spec_NEON_div_sel) <- c(
  "Site", "plotID", "MSD", "MSD2", "MSm",
  "M_s", "mPrime_s", "qHt_s", "qEt_s", "qDT_s", "qDTM_s"
)

head(spec_NEON_div_sel)

### Trait diversity
load("Results/traitDiversity/trait_DIS_tables.RData")

trait_NEON_dis_table_q0 <- trait_NEON_dis_table_q0 %>%
  select(Site, plotID, everything())

trait_NEON_dis_table_q1 <- trait_NEON_dis_table_q1 %>%
  select(Site, plotID, everything())

trait_NEON_dis_table_q2 <- trait_NEON_dis_table_q2 %>%
  select(Site, plotID, everything())

trait_NEON_dis_table_q3 <- trait_NEON_dis_table_q3 %>%
  select(Site, plotID, everything())

### Q 0
trait_NEON_dis_table_q0_sel <- trait_NEON_dis_table_q0 %>%
  select(Site, plotID, MFD, MFDz, SR, M, mPrime, qHt, qEt, qDT, qDTM)

names(trait_NEON_dis_table_q0_sel) <- c(
  "Site", "plotID", "MFD", "MFDz", "SR", "M_f",
  "mPrime_f", "qHt_f", "qEt_f", "qDT_f", "qDTM_f"
)

## combine spectral and trait diversity metrics
trait_spec_NEON_table_q0 <- full_join(trait_NEON_dis_table_q0_sel,
  spec_NEON_div_sel,
  by = c("Site", "plotID")
)

## Q 1
trait_NEON_dis_table_q1_sel <- trait_NEON_dis_table_q1 %>%
  select(Site, plotID, MFD, MFDz, SR, M, mPrime, qHt, qEt, qDT, qDTM)

names(trait_NEON_dis_table_q1_sel) <- c(
  "Site", "plotID", "MFD", "MFDz", "SR", "M_f",
  "mPrime_f", "qHt_f", "qEt_f", "qDT_f", "qDTM_f"
)

## combine spectral and trait diversity metrics
trait_spec_NEON_table_q1 <- full_join(trait_NEON_dis_table_q1_sel,
  spec_NEON_div_sel,
  by = c("Site", "plotID")
)

## Q 2
trait_NEON_dis_table_q2_sel <- trait_NEON_dis_table_q2 %>%
  select(Site, plotID, MFD, MFDz, SR, M, mPrime, qHt, qEt, qDT, qDTM)

names(trait_NEON_dis_table_q2_sel) <- c(
  "Site", "plotID", "MFD", "MFDz", "SR", "M_f",
  "mPrime_f", "qHt_f", "qEt_f", "qDT_f", "qDTM_f"
)

## combine spectral and trait diversity metrics
trait_spec_NEON_table_q2 <- full_join(trait_NEON_dis_table_q2_sel,
  spec_NEON_div_sel,
  by = c("Site", "plotID")
)

## Q 3
trait_NEON_dis_table_q3_sel <- trait_NEON_dis_table_q3 %>%
  select(Site, plotID, MFD, MFDz, SR, M, mPrime, qHt, qEt, qDT, qDTM)

names(trait_NEON_dis_table_q3_sel) <- c(
  "Site", "plotID", "MFD", "MFDz", "SR", "M_f",
  "mPrime_f", "qHt_f", "qEt_f", "qDT_f", "qDTM_f"
)

## combine spectral and trait diversity metrics
trait_spec_NEON_table_q3 <- full_join(trait_NEON_dis_table_q3_sel,
  spec_NEON_div_sel,
  by = c("Site", "plotID")
)

##### Save matched tables #####

save(trait_spec_NEON_table_q0, trait_spec_NEON_table_q1,
  trait_spec_NEON_table_q2, trait_spec_NEON_table_q3,
  file = "Results/RegDATA/trait_spec_DIS_NEON.RData"
)


##### --------------------- Trait diversity METRICS -------------------- #####

### Spectral diversity
spec_NEON_div <- readRDS("Results/spectralDiversity/spectral_NEON_q0.rds")

spec_NEON_div <- spec_NEON_div %>%
  select(Site, plotID, everything())

names(spec_NEON_div)

##### Match tables and prepare regressions #####

##### Spectral diversity
spec_NEON_div_sel <- spec_NEON_div %>%
  select(
    Site, plotID, SD, MSD, M, mPrime, qHt, qEt, qDT, qDTM, MSm, Range,
    SDivergence, SEvenness, SDispersion, SRao
  ) %>%
  mutate(SD2 = SD, MSD2 = MSD) %>%
  select(Site, plotID, MSm, SDivergence, SEvenness, SDispersion, SRao)

names(spec_NEON_div_sel) <- c(
  "Site", "plotID", "MSm", "Divergence_s",
  "Evenness_s", "Dispersion_s", "Rao_s"
)

head(spec_NEON_div_sel)

### Trait diversity
load("Results/traitDiversity/trait_DIS_tables.RData")

trait_NEON_dis_table_q0 <- trait_NEON_dis_table_q0 %>%
  select(Site, plotID, everything())

trait_NEON_dis_table_q1 <- trait_NEON_dis_table_q1 %>%
  select(Site, plotID, everything())

trait_NEON_dis_table_q2 <- trait_NEON_dis_table_q2 %>%
  select(Site, plotID, everything())

trait_NEON_dis_table_q3 <- trait_NEON_dis_table_q3 %>%
  select(Site, plotID, everything())


### Q 0
trait_NEON_table_q0_sel <- trait_NEON_dis_table_q0 %>%
  select(Site, plotID, TRichness, TDivergence, TEveness, TDispersion, TRao)

names(trait_NEON_table_q0_sel) <- c(
  "Site", "plotID", "Richness_t", "Divergence_t",
  "Evenness_t", "Dispersion_t", "Rao_t"
)

## combine spectral and trait diversity metrics
trait_spec_NEON_table_MET_q0 <- full_join(trait_NEON_table_q0_sel,
  spec_NEON_div_sel,
  by = c("Site", "plotID")
)

### Q 1
trait_NEON_table_q1_sel <- trait_NEON_dis_table_q1 %>%
  select(Site, plotID, TRichness, TDivergence, TEveness, TDispersion, TRao)

names(trait_NEON_table_q1_sel) <- c(
  "Site", "plotID", "Richness_t", "Divergence_t",
  "Evenness_t", "Dispersion_t", "Rao_t"
)

## combine spectral and trait diversity metrics
trait_spec_NEON_table_MET_q1 <- full_join(trait_NEON_table_q1_sel,
  spec_NEON_div_sel,
  by = c("Site", "plotID")
)

### Q 2
trait_NEON_table_q2_sel <- trait_NEON_dis_table_q2 %>%
  select(Site, plotID, TRichness, TDivergence, TEveness, TDispersion, TRao)

names(trait_NEON_table_q2_sel) <- c(
  "Site", "plotID", "Richness_t", "Divergence_t",
  "Evenness_t", "Dispersion_t", "Rao_t"
)

## combine spectral and trait diversity metrics
trait_spec_NEON_table_MET_q2 <- full_join(trait_NEON_table_q2_sel,
  spec_NEON_div_sel,
  by = c("Site", "plotID")
)

### Q 3
trait_NEON_table_q3_sel <- trait_NEON_dis_table_q3 %>%
  select(Site, plotID, TRichness, TDivergence, TEveness, TDispersion, TRao)

names(trait_NEON_table_q3_sel) <- c(
  "Site", "plotID", "Richness_t", "Divergence_t",
  "Evenness_t", "Dispersion_t", "Rao_t"
)

## combine spectral and trait diversity metrics
trait_spec_NEON_table_MET_q3 <- full_join(trait_NEON_table_q3_sel,
  spec_NEON_div_sel,
  by = c("Site", "plotID")
)

##### Save matched tables #####

save(trait_spec_NEON_table_MET_q0, trait_spec_NEON_table_MET_q1,
  trait_spec_NEON_table_MET_q2, trait_spec_NEON_table_MET_q3,
  file = "Results/RegDATA/trait_spec_MET_NEON.RData"
)

########### --------------- Taxonomic dimension ---------------- ############

library(tidyverse)

##### Alpha diversity #####

### Spectral species diversity
load("Results/taxoDiversity/alphaSpec_diversity_NEON.RData")



### Ground diversity
load("Results/taxoDiversity/alphaGround_diversity_NEON.RData")

### Headers
namesSpec <- c(
  "SRich_s", "SRichSErr_s", "Shannon_s", "ShannonSErr_s",
  "Simpson_s", "SimpsonSErr_s", "H_s", "Simp_s", "invSimp_s",
  "unbiasSimp_s", "Fisher_alpha_s", "S_s", "J_Pielou_s",
  "plotID", "Site"
)

namesGround <- c(
  "SRich_g", "SRichSErr_g", "Shannon_g", "ShannonSErr_g",
  "Simpson_g", "SimpsonSErr_g", "H_g", "Simp_g", "invSimp_g",
  "unbiasSimp_g", "Fisher_alpha_g", "S_g", "J_Pielou_g",
  "plotID", "Site"
)

### Rename columns
names(alphaSpec_obs_NEON_table) <- namesSpec
head(alphaSpec_obs_NEON_table)

names(alphaSpec_thresh_NEON_table) <- namesSpec
head(alphaSpec_thresh_NEON_table)

names(alphaGround_NEON_table) <- namesGround
head(alphaGround_NEON_table)

### Join diversity metrics from ground and spectra threshold
taxo_spec_thresh_NEON_table <- full_join(alphaGround_NEON_table,
  alphaSpec_thresh_NEON_table,
  by = c("Site", "plotID")
) %>%
  select(Site, plotID, everything())

### Join diversity metrics from ground and spectra threshold
taxo_spec_obs_NEON_table <- full_join(alphaGround_NEON_table,
  alphaSpec_obs_NEON_table,
  by = c("Site", "plotID")
) %>%
  select(Site, plotID, everything())

### save matched tables
save(taxo_spec_obs_NEON_table,
  taxo_spec_thresh_NEON_table,
  file = "Results/taxoDiversity/taxo_spec_ALPHA_NEON.RData"
)

##### Beta diversity #####
# Diversity on the ground
load("Results/taxoDiversity/betaGround_diversity_NEON.RData")

# Diversity estimated using spectra
load("Results/taxoDiversity/betaSpec_diversity_NEON.RData")

### Headers
namesSpec <- c("D_s", "StdErr_s", "q", "Diversity", "Site")

namesGround <- c("D_g", "StdErr_g", "q", "Diversity", "Site")

### Rename tables
names(betaGround_NEON_table) <- namesGround
head(betaGround_NEON_table)

names(betaSpec_obs_NEON_table) <- namesSpec
head(betaSpec_obs_NEON_table)

names(betaSpec_thresh_NEON_table) <- namesSpec
head(betaSpec_thresh_NEON_table)

### Join diversity tables
taxo_spec_COM_thresh_NEON_table <- full_join(betaGround_NEON_table,
  betaSpec_thresh_NEON_table,
  by = c("q", "Diversity", "Site")
) %>%
  select(Site, Diversity, q, everything())

head(taxo_spec_COM_thresh_NEON_table)

taxo_spec_COM_obs_NEON_table <- full_join(betaGround_NEON_table,
  betaSpec_obs_NEON_table,
  by = c("q", "Diversity", "Site")
) %>%
  select(Site, Diversity, q, everything())

### Save tables
save(taxo_spec_COM_thresh_NEON_table, 
     taxo_spec_COM_obs_NEON_table, 
     file = "Results/taxoDiversity/taxo_spec_BETA_NEON.RData")

##### Combine metrics based SAM #####
library(tidyverse)

##### Load and arrange data #####
### Phylo diversity
load("Results/phyloDiversity/phylo_DIS_tables.RData")

phylo_NEON_table_q0 <- phylo_NEON_table_q0 %>%
  dplyr::select(Site, plotID, everything())

phylo_NEON_table_q1 <- phylo_NEON_table_q1 %>%
  dplyr::select(Site, plotID, everything())

phylo_NEON_table_q2 <- phylo_NEON_table_q2 %>%
  dplyr::select(Site, plotID, everything())

phylo_NEON_table_q3 <- phylo_NEON_table_q3 %>%
  dplyr::select(Site, plotID, everything())

### Spectral diversity
spec_NEON_div <- readRDS("Results/spectralDiversity/spectral_NEON_SAM_q0.rds")

spec_NEON_div <- spec_NEON_div %>%
  dplyr::select(Site, plotID, everything())

##### Match tables and prepare regressions #####

##### Spectral diversity
spec_NEON_div_sel <- spec_NEON_div %>%
  dplyr::select(Site, plotID, SD, MSD, M, mPrime, qHt, qEt, qDT, qDTM, MSm, Range) %>%
  mutate(SD2 = SD, MSD2 = MSD) %>%
  dplyr::select(Site, plotID, SD, SD2, MSD, MSD2, everything())

names(spec_NEON_div_sel) <- c(
  "Site", "plotID", "SD", "SD2", "MSD", "MSD2",
  "M_s", "mPrime_s", "qHt_s", "qEt_s", "qDT_s", "qDTM_s",
  "MSm", "Range"
)

##### Phylogenetic diversity
### Q 0
names(phylo_NEON_table_q0)
names(spec_NEON_div)

## Phylogenetic diversity
phylo_NEON_table_q0_sel <- phylo_NEON_table_q0 %>%
  dplyr::select(Site, plotID, PD, PDz, MPD, MPDz, M, mPrime, qHt, qEt, qDT, qDTM)

names(phylo_NEON_table_q0_sel) <- c(
  "Site", "plotID", "PD", "PDz", "MPD", "MPDz",
  "M_p", "mPrime_p", "qHt_p", "qEt_p", "qDT_p", "qDTM_p"
)

## combine spectral and phylogenetic diversity metrics
phylo_spec_NEON_table_SAM_q0 <- full_join(phylo_NEON_table_q0_sel, spec_NEON_div_sel,
                                      by = c("Site", "plotID")
)

##### Save matched tables #####

save(phylo_spec_NEON_table_SAM_q0,
     file = "Results/RegDATA/phylo_spec_DIS_SAM_NEON.RData"
)
