library(tidyverse)
library(Taxonstand)

X <- readRDS("DATA/Traits/Traits_BIEN/traits_BIEN_NEON_spp_genus_combined.rds")

trt_world <- read.csv("../Lucie_Jesus/RCN_evolutionary_legacy/Data/TryData/FinalDATA/traitsMATCH_checkNOnas_FINAL.csv")

lpm <- trt_world %>% 
  select(Species, lpm) %>% 
  drop_na() %>% 
  group_by(Species) %>% 
  summarize(mean_lmp = mean(lpm))

sla <- trt_world %>% 
  select(Species, sla) %>% 
  drop_na() %>% 
  group_by(Species) %>% 
  summarize(mean_sla = mean(sla))

lnm <- trt_world %>% 
  select(Species, lnm) %>% 
  drop_na() %>% 
  group_by(Species) %>% 
  summarize(mean_lnm = mean(lnm))

lpm <- data.frame(lpm)
sla <- data.frame(sla)
lnm <- data.frame(lnm)

traits <- full_join(X, lpm, by = c("Taxa" = "Species"))
traits <- full_join(traits, sla, by = c("Taxa" = "Species"))
traits <- full_join(traits, lnm, by = c("Taxa" = "Species"))

### Select traits 
sels <- seq(1, 33, by = 2)

trait_sel <- traits[4:36][, sels]
trait_sel <- cbind(traits[1:3], trait_sel)
trait_sel <- cbind(trait_sel, traits[37:38])
## Taxonomy 
spp <- unique(trait_sel$Taxa) # vector with scientific names
spp <- gsub("_", " ", spp)
# Perform taxonomic standardization on plant names (TPL table)
spp_check <- TPL(spp, infra = TRUE, corr = TRUE)

# Select the necessary information and combine with the trait data 
load("DATA/Traits/Traits_BIEN/GapFilling/taxonomy_traits_update.RData")

taxonomy <- spp_check %>% 
  mutate("species2" = paste(New.Genus, New.Species)) %>%  
  drop_na(New.Genus, Taxon, Family) 

taxonomy <- taxonomy %>% 
  mutate("plant_id" = 1:nrow(taxonomy), "species" = Taxon, 
         "genus" = New.Genus, "family" = Family) %>% 
  select(plant_id, species, genus, family)

tail(taxonomy) 

gen_fam <- read.csv("DATA/Traits/Traits_BIEN/GapFilling/taxonomy/gen_family_clean.csv")

taxonomy2 <- full_join(taxonomy, gen_fam[, 2:3], 
                       by = c("genus" = "genus"))

taxonomy2 <- taxonomy2 %>% 
  distinct(species, .keep_all = TRUE)

sum(complete.cases(taxonomy2))
sum(!complete.cases(taxonomy2))

#save(spp_check, taxonomy, 
 #    file = "DATA/Traits/Traits_BIEN/GapFilling/taxonomy_traits_update.RData")

trait_sel <- trait_sel[, c(1, 4:22)] 

trait_sel <- trait_sel %>% 
  mutate(species = gsub("_", " ", Taxa)) %>% 
  select(species, everything()) 

head(trait_sel)

trait_sel <- trait_sel[, c(1, 3:21)]

trait_imputation <- full_join(taxonomy2, trait_sel, 
                              by = "species" )

apply(!is.na(trait_imputation), 2, sum)


trait_imputation$family[trait_imputation$family == ""] <- NA

trait_imputation <- trait_imputation %>% 
  drop_na(family.x)

trait_imputation$family[trait_imputation$family == ""]
trait_imputation$genus[trait_imputation$genus == ""]
trait_imputation$species[trait_imputation$species == ""]

sum(complete.cases(trait_imputation[1:5]))
sum(!complete.cases(trait_imputation[1:5]))


trait_imputation <- trait_imputation[order(trait_imputation$family.x, 
                                           trait_imputation$genus), ] 

apply(!is.na(trait_imputation), 2, sum)

#write.csv(trait_imputation, 
 #         file = "DATA/Traits/Traits_BIEN/GapFilling/trait_to_imputation.csv")

write.csv(trait_imputation, 
          file = "DATA/Traits/Traits_BIEN/GapFilling/trait_to_imputation_update.csv") 

phyNEON <- read.nexus("DATA/Phylogeny/Plants/phylo_ALL_NEON_S3.nex")

trait_imputation$species <- gsub(" ", "_", trait_imputation$species)

traitPhylo <- trait_imputation[trait_imputation$species %in% phyNEON$tip.label, ]

apply(!is.na(traitPhylo), 2, sum)


##### Prepare data for imputation #####
library(tidyverse)

trait_check <- read.csv("DATA/Traits/Traits_BIEN/GapFilling/trait_to_imputation_update.csv")

trait_check2 <- trait_check %>% 
  mutate(SR = rowSums(trait_check[, 7:25], na.rm = TRUE) * NA ^ (rowSums(!is.na(trait_check[, 7:25])) == 0)
  ) %>% 
  drop_na(SR)

trait_check2$plant_id <- 1:nrow(trait_check2)

hierarchy.info <- trait_check2[, c(2:4, 6)]

hierarchy.info$species <- as.factor(hierarchy.info$species)
hierarchy.info$genus <- as.factor(hierarchy.info$genus)
hierarchy.info$family <- as.factor(hierarchy.info$family.x)

hierarchy.info <- hierarchy.info %>% 
  select(plant_id, species, genus, family)

head(hierarchy.info, 10)
tail(hierarchy.info, 10)

str(hierarchy.info)

apply(is.na(hierarchy.info), 2, sum)
apply(!is.na(hierarchy.info), 2, sum)

## Check for NA
hierarchy.info %>% 
  filter(is.na(plant_id))

hierarchy.info %>% 
  filter(is.na(species))

hierarchy.info %>% 
  filter(is.na(genus))

hierarchy.info %>% 
  filter(is.na(family))

## Trait values
trait.info <- trait_check2[, c(7:25)]
#names(trait.info) <- c("Trait1", "Trait2", "Trait3", "Trait4", 
 #                      "Trait5", "Trait6", "Trait7", "Trait8", 
  #                     "Trait9", "Trait10", "Trait11", "Trait12", 
   #                    "Trait13", "Trait14", "Trait15", "Trait16", "Trait17")

trait.info <- as.matrix(trait.info)
dim(trait.info)
head(trait.info, 10)
tail(trait.info, 10)

#sum(row.names(trait.info) %in% colnames(hierarchy.info)) == ncol(hierarchy.info)
nrow(hierarchy.info) == nrow(trait.info)

trait.info2 <- apply(X = trait.info, MARGIN = 2, FUN = scale)

trait.info3 <- apply(X = trait.info, MARGIN = 2, FUN = log10)

nrow(hierarchy.info) == nrow(trait.info2)

sum(!is.na(trait.info[, 17]))
sum(is.na(trait.info[, 17]))

apply(is.na(trait.info), 2, sum)
apply(!is.na(trait.info), 2, sum)

apply(!is.na(trait.info), 1, sum)

##### run Bayesian gapfilling #####
library(BHPMF)

GapFilling(X = trait.info2, 
           hierarchy.info = hierarchy.info, 
           num.samples = 10000, burn = 2000, 
           tuning = TRUE, verbose = TRUE, 
           tmp.dir = "DATA/Traits/Traits_BIEN/GapFilling/tmp7", 
           mean.gap.filled.output.path = "DATA/Traits/Traits_BIEN/GapFilling/tmp7/NEON_mean_gap_filled.csv",
           std.gap.filled.output.path = "DATA/Traits/Traits_BIEN/GapFilling/tmp7/NEON_std_gap_filled.csv")

#RRMSE for the test data:  1.205673$min.rmse
#[1] 1.370784

#$best.number.latent.features
#[1] 20

##########################################
# Usage 4: Calculate cross validation RMSE
##########################################
# Calculate average RMSE with the default values 
traitNames <- colnames(trait.info2)

for(i in 1:ncol(trait.info2)) { 
  
  print(traitNames[i])
  
  dir.create(paste0("DATA/Traits/Traits_BIEN/GapFilling/CV/", traitNames[i], "_cv"))
  
  out1 <- CalculateCvRmse(as.matrix(trait.info2[, i]), 
                          hierarchy.info, 
                          #tuning = TRUE, 
                          num.latent.feats = 20, 
                          tmp.dir = paste0("DATA/Traits/Traits_BIEN/GapFilling/CV/", traitNames[i], "_cv"), 
                          num.samples = 1000, 
                          burn = 200, 
                          verbose = TRUE)
  
  #avg.rmse <- out1$avg.rmse
  #std.rmse <- out1$std.rmse
  
  save(out1, file = paste0("DATA/Traits/Traits_BIEN/GapFilling/CV/", traitNames[i], "_cv.RData"))
  
  
  }

out1 <- CalculateCvRmse(trait.info2, 
                        hierarchy.info, 
                        #tuning = TRUE, 
                        num.latent.feats = 20, 
                        tmp.dir = "DATA/Traits/Traits_BIEN/GapFilling/tmp_cv6", 
                        num.samples = 10000, 
                        burn = 2000, 
                        verbose = TRUE)

avg.rmse <- out1$avg.rmse
std.rmse <- out1$std.rmse

save(out1, file = "DATA/Traits/Traits_BIEN/GapFilling/CV/cv7_rmse.RData")



GapFilling(X = trait.info2, 
           hierarchy.info = hierarchy.info, 
           num.samples = 10000, burn = 2000,
           tmp.dir = "DATA/Traits/Traits_BIEN/GapFilling/tmp2", 
           mean.gap.filled.output.path = "DATA/Traits/Traits_BIEN/GapFilling/NEON_mean_gap_filled_noScale.csv",
           std.gap.filled.output.path = "DATA/Traits/Traits_BIEN/GapFilling/NEON_std_gap_filled_noScale.csv")

out2 <- CalculateCvRmse(trait.info, 
                        hierarchy.info, 
                        num.samples = 10000, 
                        burn = 2000)

save(out2, file = "DATA/Traits/Traits_BIEN/GapFilling/cv_rmse_noScale.RData")

## Gapfilling using log

GapFilling(X = trait.info3, 
           hierarchy.info = hierarchy.info, 
           num.samples = 10000, burn = 2000,
           tmp.dir = "DATA/Traits/Traits_BIEN/GapFilling/tmp3", 
           mean.gap.filled.output.path = "DATA/Traits/Traits_BIEN/GapFilling/NEON_mean_gap_filled_log.csv",
           std.gap.filled.output.path = "DATA/Traits/Traits_BIEN/GapFilling/NEON_std_gap_filled_log.csv",
           verbose = TRUE)

##########################################
# Usage 4: Calculate cross validation RMSE
##########################################
# Calculate average RMSE with the default values
out3 <- CalculateCvRmse(trait.info3, 
                        hierarchy.info, 
                        tuning = TRUE, 
                        tmp.dir = "DATA/Traits/Traits_BIEN/GapFilling/tmp_cv3", 
                        num.samples = 10000, 
                        burn = 2000, 
                        verbose = TRUE)

out3

save(out3, file = "DATA/Traits/Traits_BIEN/GapFilling/cv_rmse_log.RData")


## Not run: 
# Read the input data
data(hierarchy.info) # Read the hierarchy information
head(hierarchy.info)
str(hierarchy.info)

data(trait.info)  # Read the matrix X
head(trait.info)
str(trait.info)

trait.info2 <- apply(trait.info, 2, scale)
  
sum(row.names(trait.info) %in% colnames(hierarchy.info)) == ncol(hierarchy.info)
nrow(hierarchy.info) == nrow(trait.info)

back_trans_pars <- list()
rm_col <- c()
for(i in 1:ncol(trait.info)){
  x <- trait.info[,i] # goes through the columns
  min_x <- min(x,na.rm = T) # takes the min of each column
  if(min_x < 0.00000000001){
    x <- x - min_x + 1 # make this optional if min x is neg
  }
  logx <- log10(x)
  mlogx <- mean(logx, na.rm = T)
  slogx <- sd(logx, na.rm = T)
  x <- (logx - mlogx)/slogx # Z transformation
  back_trans_pars[[i]] <- list(min_x = min_x,
                               mlogx = mlogx,
                               slogx = slogx)
  trait.info[,i] <- x
}
write.csv(back_trans_pars, paste0("Downloads/back_trans_pars.csv"))


###################
# Usage 1: Figure 1
###################
GapFilling(trait.info, hierarchy.info,
           mean.gap.filled.output.path = "../../Downloads/mean_gap_filled.txt",
           std.gap.filled.output.path="../../Downloads/std_gap_filled.txt")

###################
# Usage 2: Figure 2
###################
GapFilling(trait.info, hierarchy.info,
           prediction.level = 4,
           used.num.hierarchy.levels = 2,
           mean.gap.filled.output.path = "../../Downloads/mean_gap_filled.txt",
           std.gap.filled.output.path="/../../Downloads/std_gap_filled.txt")

###################
# Usage 3: Figure 3
###################
GapFilling(trait.info, hierarchy.info,
           prediction.level = 3,
           used.num.hierarchy.levels = 2,
           mean.gap.filled.output.path = "../../Downloads/mean_gap_filled.txt",
           std.gap.filled.output.path="../../Downloads/std_gap_filled.txt")  

######################
# Usage 4: Running PMF
# not using hierarchy
#####################
GapFilling(trait.info, hierarchy.info,
           used.num.hierarchy.levels = 0,
           mean.gap.filled.output.path = "../../Downloads/mean_gap_filled.txt",
           std.gap.filled.output.path="../../Downloads/std_gap_filled.txt")  


##########################################
# Usage 4: Calculate cross validation RMSE
##########################################
# Calculate average RMSE with the default values
out1 <- CalculateCvRmse(trait.info, hierarchy.info)
avg.rmse <- out1$avg.rmse
std.rmse <- out1$std.rmse

# Calculate average RMSE using partial hierarch: using only the first 2 level of hierarchy
out2 <- CalculateCvRmse(X, hierarchy.info, used.num.hiearchy.levels=2,
                        num.samples=500, burn=20, gaps=3, tuning=TRUE)
avg.rmse <- out2$avg.rmse
std.rmse <-	out2$std.rmse

## End(Not run)