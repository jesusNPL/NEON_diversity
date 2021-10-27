##### Bayesian correlation taxonomic dimension #####

demon_BayCorrTAX <- function(matRES, nChains, nIters, nCores, pathSave) {
  library(brms)
  library(rstan)
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  # Set prior for Bayesian robust correlations
  Corr_prior <- c(prior(gamma(2, 0.1), class = nu),
                  prior(lkj(1), class = rescor))
  
  fits <- list()
  R2_robust <- list()
  R2_regular <- list()
  
  for(i in 1:10) {
    print(names(matRES[i + 1]))
    print(names(matRES[i + 11]))
    
    data <- matRES[, c(i+1, i+11)] 
    
    headers <- names(data)
    
    formulas <- as.formula(paste(headers[1], " ~ ", paste(headers[2], collapse = "+"))) 
    
    fit <- brm(formulas, data = data, 
               #family = student(), 
               iter = nIters, 
               warmup = nIters/5, 
               chains = nChains, 
               cores = nCores) 
    robust <- data.frame(bayes_R2(fit, robust = TRUE))
    robust$Veg <- headers[1]
    robust$Spec <- headers[2]
    R2_robust[[i]] <- robust
    
    regular <- data.frame(bayes_R2(fit))
    regular$Veg <- headers[1]
    regular$Spec <- headers[2]
    R2_regular[[i]] <- regular
    
    fits[[i]] <- fit 
  }
  
  R2_robust <- do.call(rbind, R2_robust)
  R2_regular <- do.call(rbind, R2_regular)
  
  results <- list(CorrRobust = R2_robust, Corr = R2_regular, fits = fits)
  
  save(results, file = pathSave)
  
  return(results)
  
}

##### Bayesian correlation phylogenetic dimension #####

demon_BayCorrPHY <- function(matRES, nChains, nIters, nCores, pathSave, Q) {
  library(brms)
  library(rstan)
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  # Set prior for Bayesian robust correlations
  Corr_prior <- c(prior(gamma(2, 0.1), class = nu),
                  prior(lkj(1), class = rescor)) 
  mat <- matRES %>% 
    drop_na()
  
  mat2 <- mat %>% 
    filter(qHill == Q) %>% 
    mutate(SDz = SD, MSDz = MSD) %>% 
    select(plotID, PD, PDz, MPD, MPDz, Rao_Phylo, 
           M_Phylo, mPrime_Phylo, qHt_Phylo, qEt_Phylo, qDT_Phylo, qDTM_Phylo, 
           SD, SDz, MSD, MSDz, Rao_Spec, 
           M_Spec, mPrime_Spec, qHt_Spec, qEt_Spec, qDT_Spec, qDTM_Spec)
  
  fits <- list()
  R2_robust <- list()
  R2_regular <- list()
  
  for(i in 1:11) {
    print(names(mat2[i + 1]))
    print(names(mat2[i + 12]))
    
    data <- mat2[, c(i + 1, i + 12)] 
    
    headers <- names(data)
    
    formulas <- as.formula(paste(headers[1], " ~ ", paste(headers[2], collapse = "+"))) 
    
    fit <- brm(formulas, data = data, 
               #family = student(), 
               iter = nIters, 
               warmup = nIters/5, 
               chains = nChains, 
               cores = nCores) 
    robust <- data.frame(bayes_R2(fit, robust = TRUE))
    robust$Phylo <- headers[1]
    robust$Spec <- headers[2]
    R2_robust[[i]] <- robust
    
    regular <- data.frame(bayes_R2(fit))
    regular$Phylo <- headers[1]
    regular$Spec <- headers[2]
    R2_regular[[i]] <- regular
    
    fits[[i]] <- fit 
  }
  
  R2_robust <- do.call(rbind, R2_robust)
  R2_regular <- do.call(rbind, R2_regular)
  
  results <- list(CorrRobust = R2_robust, Corr = R2_regular, fits = fits)
  
  save(results, file = pathSave)
  
  return(results)
  
  print(paste0("Bayesian models for Q = ", Q, " done..."))
  
}

##### Bayesian correlation trait dimension #####

demon_BayCorrTRAIT <- function(matRES, nChains, nIters, nCores, pathSave, Q) {
  library(brms)
  library(rstan)
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  # Set prior for Bayesian robust correlations
  Corr_prior <- c(prior(gamma(2, 0.1), class = nu),
                  prior(lkj(1), class = rescor)) 
  
  mat2 <- matRES %>% 
    filter(qHill == Q) %>% 
    mutate(MSDz = MSD) %>% 
    select(plotID, MFD, MFDz, Rao_Trait, 
           M_Trait, mPrime_Trait, qHt_Trait, qEt_Trait, qDT_Trait, qDTM_Trait, 
           Divergence_Trait, Evenness_Trait, Dispersion_Trait, 
           MSD, MSDz, Rao_Spec, 
           M_Spec, mPrime_Spec, qHt_Spec, qEt_Spec, qDT_Spec, qDTM_Spec, 
           Divergence_Spec, Evenness_Spec, Dispersion_Spec)
  
  fits <- list()
  R2_robust <- list()
  R2_regular <- list()
  
  for(i in 1:12) {
    print(names(mat2[i + 1]))
    print(names(mat2[i + 13]))
    
    data <- mat2[, c(i + 1, i + 13)] 
    
    headers <- names(data)
    
    formulas <- as.formula(paste(headers[1], " ~ ", paste(headers[2], collapse = "+"))) 
    
    fit <- brm(formulas, data = data, 
               #family = student(), 
               iter = nIters, 
               warmup = nIters/5, 
               chains = nChains, 
               cores = nCores) 
    robust <- data.frame(bayes_R2(fit, robust = TRUE))
    robust$Trait <- headers[1]
    robust$Spec <- headers[2]
    R2_robust[[i]] <- robust
    
    regular <- data.frame(bayes_R2(fit))
    regular$Trait <- headers[1]
    regular$Spec <- headers[2]
    R2_regular[[i]] <- regular
    
    fits[[i]] <- fit 
  }
  
  R2_robust <- do.call(rbind, R2_robust)
  R2_regular <- do.call(rbind, R2_regular)
  
  results <- list(CorrRobust = R2_robust, Corr = R2_regular, fits = fits)
  
  save(results, file = pathSave)
  
  return(results)
  
  print(paste0("Bayesian models for Q = ", Q, " done..."))
}

###### Bayesian moments #####
demon_BayCorrMOMENTS <- function(momentsRES, nChains, nIters, nCores, 
                                 pathSaveTaxo, pathSaveTrait, pathSavePhylo) {
  library(brms)
  library(rstan)
  
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  # Set prior for Bayesian robust correlations
  Corr_prior <- c(prior(gamma(2, 0.1), class = nu),
                  prior(lkj(1), class = rescor)) 
  
  taxo <- momentsRES$taxoSpec
  trait <- momentsRES$traitSpec
  phylo <- momentsRES$phyloSpec

  # Taxonomy
  fits_taxo <- list()
  R2_robust_taxo <- list()
  R2_regular_taxo <- list()
  # Trait
  fits_trait <- list()
  R2_robust_trait <- list()
  R2_regular_trait <- list()
  # Phylogeny
  fits_phylo <- list()
  R2_robust_phylo <- list()
  R2_regular_phylo <- list()
  
  for(i in 1:6) {
    print("Taxonomy ~ Spectra")
    print(names(taxo[i + 1]))
    print(names(taxo[i + 7]))
    
    dataTx <- taxo[, c(i + 1, i + 7)] 
    headersTx <- names(dataTx)
    formTx <- as.formula(paste(headersTx[1], " ~ ", 
                               paste(headersTx[2], collapse = "+"))) 
    
    fitTx <- brm(formTx, data = dataTx, 
               #family = student(), 
               iter = nIters, 
               warmup = nIters/5, 
               chains = nChains, 
               cores = nCores) 
    
    robustTaxo <- data.frame(bayes_R2(fitTx, robust = TRUE))
    robustTaxo$Taxo <- headersTx[1]
    robustTaxo$Spec <- headersTx[2]
    R2_robust_taxo[[i]] <- robustTaxo
    
    regularTaxo <- data.frame(bayes_R2(fitTx))
    regularTaxo$Taxo <- headersTx[1]
    regularTaxo$Spec <- headersTx[2]
    R2_regular_taxo[[i]] <- regularTaxo
    
    fits_taxo[[i]] <- fitTx 
    
    print("Trait ~ Spectra")
    print(names(trait[i + 1]))
    print(names(trait[i + 7]))
    
    dataTrt <- trait[, c(i + 1, i + 7)] 
    headersTrt <- names(dataTrt)
    formTrt <- as.formula(paste(headersTrt[1], " ~ ", paste(headersTrt[2], collapse = "+"))) 
    
    fitTrt <- brm(formTrt, data = dataTrt, 
                 #family = student(), 
                 iter = nIters, 
                 warmup = nIters/5, 
                 chains = nChains, 
                 cores = nCores) 
    
    robustTrait <- data.frame(bayes_R2(fitTrt, robust = TRUE))
    robustTrait$Trait <- headersTrt[1]
    robustTrait$Spec <- headersTrt[2]
    R2_robust_trait[[i]] <- robustTrait
    
    regularTrait <- data.frame(bayes_R2(fitTrt))
    regularTrait$Taxo <- headersTrt[1]
    regularTrait$Spec <- headersTrt[2]
    R2_regular_trait[[i]] <- regularTrait
    
    fits_trait[[i]] <- fitTrt 
    
    print("Phylogeny ~ Spectra")
    print(names(phylo[i + 1]))
    print(names(phylo[i + 7]))
    
    dataPhy <- phylo[, c(i + 1, i + 7)] 
    headersPhy <- names(dataPhy)
    formPhy <- as.formula(paste(headersPhy[1], " ~ ", paste(headersPhy[2], collapse = "+"))) 
    
    fitPhy <- brm(formPhy, data = dataPhy, 
                  #family = student(), 
                  iter = nIters, 
                  warmup = nIters/5, 
                  chains = nChains, 
                  cores = nCores) 
    
    robustPhylo <- data.frame(bayes_R2(fitPhy, robust = TRUE))
    robustPhylo$Phylo <- headersPhy[1]
    robustPhylo$Spec <- headersPhy[2]
    R2_robust_phylo[[i]] <- robustPhylo
    
    regularPhylo <- data.frame(bayes_R2(fitPhy))
    regularPhylo$Phylo <- headersPhy[1]
    regularPhylo$Spec <- headersPhy[2]
    R2_regular_phylo[[i]] <- regularPhylo
    
    fits_phylo[[i]] <- fitPhy 
    
    
  }
  # Taxonomy
  R2_robust_taxo <- do.call(rbind, R2_robust_taxo)
  R2_regular_taxo <- do.call(rbind, R2_regular_taxo)
  
  resultsTaxo <- list(CorrRobust = R2_robust_taxo, 
                      Corr = R2_regular_taxo, 
                      fits = fits_taxo)
  save(resultsTaxo, file = pathSaveTaxo)
  # Trait
  R2_robust_trait <- do.call(rbind, R2_robust_trait)
  R2_regular_trait <- do.call(rbind, R2_regular_trait)
  
  resultsTrait <- list(CorrRobust = R2_robust_trait, 
                       Corr = R2_regular_trait, 
                       fits = fits_trait)
  save(resultsTrait, file = pathSaveTrait)
  # Phylogeny
  R2_robust_phylo <- do.call(rbind, R2_robust_phylo)
  R2_regular_phylo <- do.call(rbind, R2_regular_phylo)
  
  resultsPhylo <- list(CorrRobust = R2_robust_phylo, 
                       Corr = R2_regular_phylo, 
                       fits = fits_phylo)
  save(resultsPhylo, file = pathSavePhylo)
  #return(results)
  
  print("Bayesian models for multiple moments... done!!!")
}
