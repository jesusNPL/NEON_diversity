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

