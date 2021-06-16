##### Function to estimate spectral moments #####
## spectra = spectra data
## nPlots = number of plots
## plotNames = names of the plots
## specRange = position of the bands in the data frame.

demon_momentsSPEC <- function(spectra, nPlots, plotNames, specRange) {
  require(moments)
  library(dplyr)
  library(tidyr)
  
  ## Store calculations 
  meanSpec <- numeric(length = nPlots)
  medianSpec <- numeric(length = nPlots)
  sdSpec <- numeric(length = nPlots)
  varSpec <- numeric(length = nPlots)
  
  skew_mean_Spec <- numeric(length = nPlots)
  kurt_mean_Spec <- numeric(length = nPlots)
  skew_median_Spec <- numeric(length = nPlots)
  kurt_median_Spec <- numeric(length = nPlots)
  
  ## Start calculations
  for(i in 1:nPlots) {
    svMisc::progress(i, max.value = nPlots)
    
    spec <- spectra %>% filter(plotID == plotNames[i]) %>% 
      na.omit()

    spec <- spec[, specRange]
    
    meanSpec[i] <- mean(abs(na.omit(apply(abs(spec), MARGIN = 2, FUN = mean))))
    medianSpec[i] <- mean(abs(na.omit(apply(abs(spec), MARGIN = 2, FUN = median))))
    sdSpec[i] <- mean(abs(na.omit(apply(abs(spec), MARGIN = 2, FUN = sd))))
    varSpec[i] <- mean(abs(na.omit(apply(abs(spec), MARGIN = 2, FUN = var))))
    skew_mean_Spec[i] <- mean(abs(na.omit(apply(abs(spec), MARGIN = 2, FUN = skewness))))
    kurt_mean_Spec[i] <- mean(abs(na.omit(apply(abs(spec), MARGIN = 2, FUN = kurtosis))))
    skew_median_Spec[i] <- median(abs(na.omit(apply(abs(spec), MARGIN = 2, FUN = skewness))))
    kurt_median_Spec[i] <- median(abs(na.omit(apply(abs(spec), MARGIN = 2, FUN = kurtosis))))
    
  }
  
  momentsSpec <- data.frame(plotNames, meanSpec, medianSpec, sdSpec, varSpec, 
                            skew_mean_Spec, kurt_mean_Spec, skew_median_Spec, kurt_median_Spec) 
  
  names(momentsSpec) <- c("plotID", "mean_Spec", "median_Spec", "sd_Spec", "var_Spec", 
                          "skew_Spec", "kurt_Spec", "skew_median_Spec", "kurt_median_Spec")
  
  return(momentsSpec)
}


##### Function to estimate trait and phylogenetic moments #####
## comm = community data
## trait = trait values
## phylo = phylogenetic tree 
## nPlots = number of plots
## plotNames = names of the plots
## specRange = position of the bands in the data frame.

demon_momentsTraitPhylo <- function(comm, trait, phylo, plotNames, nPlots) {
  require(moments)
  require(picante)
  
  ## Store calculations taxonomy
  meanTaxo <- numeric(length = nPlots)
  medianTaxo <- numeric(length = nPlots)
  sdTaxo <- numeric(length = nPlots)
  varTaxo <- numeric(length = nPlots)
  
  skewTaxo <- numeric(length = nPlots)
  kurtTaxo <- numeric(length = nPlots)
  
  ## Store calculations trait
  meanTrait <- numeric(length = nPlots)
  medianTrait <- numeric(length = nPlots)
  sdTrait <- numeric(length = nPlots)
  varTrait <- numeric(length = nPlots)
  
  skewTrait <- numeric(length = nPlots)
  kurtTrait <- numeric(length = nPlots)
  
  ## Store calculations phylogeny
  phyDat <- picante::evol.distinct(tree = phylo, type = "equal.splits")
  
  meanPhylo <- numeric(length = nPlots)
  medianPhylo <- numeric(length = nPlots)
  sdPhylo <- numeric(length = nPlots)
  varPhylo <- numeric(length = nPlots)
  
  skewPhylo <- numeric(length = nPlots)
  kurtPhylo <- numeric(length = nPlots)
  
  ## Start calculations
  comm2 <- as.matrix(comm)
  
  for(i in 1:nPlots) {
    spNames <- comm2[i, ] > 0 
    selection <- names(spNames[spNames == TRUE])
    
    trait_comm <- trait[trait$Species %in% selection, ]
    phylo_comm <- phyDat[phyDat$Species %in% selection, ]
    
    ## Taxonomic moments
    comm_tmp <- t(comm[i, ])
    colnames(comm_tmp) <- "Abun"
    
    comm_tmp <- comm_tmp[rownames(comm_tmp) %in% selection, ]
    
    meanTaxo[i] <- mean(na.omit(comm_tmp))
    medianTaxo[i] <- median(na.omit(comm_tmp))
    sdTaxo[i] <- sd(na.omit(comm_tmp))
    varTaxo[i] <- var(na.omit(comm_tmp))
    
    skewTaxo[i] <- skewness(na.omit(comm_tmp))
    kurtTaxo[i] <- kurtosis(na.omit(comm_tmp))
    
    ## Trait moments
    meanTrait[i] <- mean(na.omit(trait_comm[, 2]))
    medianTrait[i] <- median(na.omit(trait_comm[, 2]))
    sdTrait[i] <- sd(na.omit(trait_comm[, 2]))
    varTrait[i] <- var(na.omit(trait_comm[, 2]))
    skewTrait[i] <- skewness(na.omit(trait_comm[, 2]))
    kurtTrait[i] <- kurtosis(na.omit(trait_comm[, 2]))
    
    ## Phylogenetic moments
    meanPhylo[i] <- mean(na.omit(phylo_comm[, 2]))
    medianPhylo[i] <- median(na.omit(phylo_comm[, 2]))
    sdPhylo[i] <- sd(na.omit(phylo_comm[, 2]))
    varPhylo[i] <- var(na.omit(phylo_comm[, 2]))
    skewPhylo[i] <- skewness(na.omit(phylo_comm[, 2]))
    kurtPhylo[i] <- kurtosis(na.omit(phylo_comm[, 2]))
  } 
  # Taxonomy
  momentsTaxo <- data.frame(plotNames, meanTaxo, medianTaxo, sdTaxo, varTaxo, 
                              skewTaxo, kurtTaxo) 
  names(momentsTaxo) <- c("plotID", "mean_Taxo", "median_Taxo", "sd_Taxo", 
                            "var_Taxo", "skew_Taxo", "kurt_Taxo")
  # Trait
  momentsTraits <- data.frame(plotNames, meanTrait, medianTrait, sdTrait, varTrait, 
                            skewTrait, kurtTrait) 
  names(momentsTraits) <- c("plotID", "mean_Trait", "median_Trait", "sd_Trait", 
                          "var_Trait", "skew_Trait", "kurt_Trait")
  # Phylogeny
  momentsPhylo <- data.frame(plotNames, meanPhylo, medianPhylo, sdPhylo, varPhylo, 
                              skewPhylo, kurtPhylo) 
  names(momentsPhylo) <- c("plotID", "mean_Phylo", "median_Phylo", "sd_Phylo", 
                            "var_Phylo", "skew_Phylo", "kurt_Phylo")
  
  momentos <- list(moments_taxo = momentsTaxo, 
                   moments_traits = momentsTraits, 
                   moments_phylo = momentsPhylo)
  return(momentos)

}


