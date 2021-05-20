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
