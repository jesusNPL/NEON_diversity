### Function to run regressions between simple distance metrics
god_BayReg_trait_dis <- function(resMetrics, 
                                 Q, 
                                 nMetrics,
                                 pathSave, 
                                 scaled = FALSE,
                                 nChains, 
                                 nIters, 
                                 burnin, 
                                 nCores, 
                                 control, 
                                 engine) { 
  
  ### Main library
  library(brms)

  ### Store results
  fits <- list()
  R2_robust <- list()
  R2_regular <- list()

  for (i in 1:nMetrics) { 
    
    print(paste0("Starting BLM for ..."))
    
    print(names(resMetrics[i + 2]))
    print(names(resMetrics[i + 6]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 6)]
    headers <- names(data)

    if (scaled == FALSE) { 
      
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste(headers[4], paste0("+ (1|siteID)"),
          collapse = "+"
        )
      )) 
      
    } else { 
      
      formulas <- as.formula(paste(headers[3],  " ~ ",
        paste0("scale(", headers[4], ")"),
        paste0("+ (1|siteID)"),
        collapse = "+"
      )) 
      
    }

    ### Run
    fit <- brms::brm(
      formula = formulas,
      data = data,
      iter = nIters,
      warmup = burnin,
      cores = nCores, 
      control = control, 
      backend = engine
    )
    
    ## extract robust Bayes R2
    robust <- data.frame(bayes_R2(fit, robust = TRUE))
    robust$Trait <- headers[3]
    robust$Spec <- headers[4] 
    robust$Q <- Q
    R2_robust[[i]] <- robust

    ## extract Bayes R2
    reg <- data.frame(bayes_R2(fit))
    reg$Trait <- headers[3]
    reg$Spec <- headers[4] 
    reg$Q <- Q
    R2_regular[[i]] <- reg

    fits[[i]] <- fit
    
  }

  ## Combine results
  R2_robust <- do.call(rbind, R2_robust)
  R2_regular <- do.call(rbind, R2_regular)

  ## Create a list to save results
  res <- list(fits = fits, R2_robust = R2_robust, R2 = R2_regular)

  save(res, file = pathSave)

  return(res)

  print(paste0("Bayesian models for Q = ", Q, " done...")) 
  
}


### Function to run regressions between classic metrics
god_BayReg_trait_met <- function(resMetrics, Q, nMetrics = 5,
                                 pathSave, scaled = FALSE,
                                 nChains, nIters, nCores, engine) {
  library(brms)
  library(cmdstanr)
  library(rstan)

  ### Inits
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  ### Store results
  fits <- list()
  R2_robust <- list()
  R2_regular <- list()

  for (i in 1:nMetrics) {
    print(names(resMetrics[i + 2]))
    print(names(resMetrics[i + 7]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 7)]
    headers <- names(data)

    if (scaled == FALSE) {
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste(headers[4], paste0("+ (1|Site)"),
          collapse = "+"
        )
      ))
    } else {
      formulas <- as.formula(paste(headers[3],  " ~ ",
        paste0("scale(", headers[4], ")"),
        paste0("+ (1|Site)"),
        collapse = "+"
      ))
    }

    ### Run
    fit <- brms::brm(
      formula = formulas,
      data = data,
      iter = nIters,
      warmup = nIters / 5,
      cores = nCores,
      backend = engine
    )
    ## extract robust Bayes R2
    robust <- data.frame(bayes_R2(fit, robust = TRUE))
    robust$Trait <- headers[3]
    robust$Spec <- headers[4]
    R2_robust[[i]] <- robust

    ## extract Bayes R2
    reg <- data.frame(bayes_R2(fit))
    reg$Trait <- headers[3]
    reg$Spec <- headers[4]
    R2_regular[[i]] <- reg

    fits[[i]] <- fit
  }

  R2_robust <- do.call(rbind, R2_robust)
  R2_regular <- do.call(rbind, R2_regular)

  res <- list(fits = fits, R2_robust = R2_robust, R2 = R2_regular)

  save(res, file = pathSave)

  return(res)

  print(paste0("Bayesian models for Q = ", Q, " done..."))
}

############ ------------- Multilevel 2 ----------------- #################

### Function to run regressions between simple distance metrics
god_BayReg_trait_dis_ML2 <- function(resMetrics, Q, nMetrics, pathSave,
                                     nChains, nIters, nCores, engine) {
  library(brms)
  library(cmdstanr)
  library(rstan)

  ### Inits
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  ### Store results
  fits <- list()
  R2_robust <- list()
  R2_regular <- list()

  for (i in 1:nMetrics) {
    print(names(resMetrics[i + 2]))
    print(names(resMetrics[i + 11]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 11)]
    headers <- names(data)

    formulas <- as.formula(paste(
      headers[3],
      " ~ ",
      paste(headers[4],
        paste0("+ (1|plotID)"),
        paste0("+ (1|Site)"),
        collapse = "+"
      )
    ))

    ### Run
    fit <- brms::brm(
      formula = formulas,
      data = data,
      iter = nIters,
      warmup = nIters / 5,
      cores = nCores,
      backend = engine
    )
    ## extract robust Bayes R2
    robust <- data.frame(bayes_R2(fit, robust = TRUE))
    robust$Trait <- headers[3]
    robust$Spec <- headers[4]
    R2_robust[[i]] <- robust

    ## extract Bayes R2
    reg <- data.frame(bayes_R2(fit))
    reg$Trait <- headers[3]
    reg$Spec <- headers[4]
    R2_regular[[i]] <- reg

    fits[[i]] <- fit
  }

  R2_robust <- do.call(rbind, R2_robust)
  R2_regular <- do.call(rbind, R2_regular)

  res <- list(fits = fits, R2_robust = R2_robust, R2 = R2_regular)

  save(res, file = pathSave)

  return(res)

  print(paste0("Bayesian models for Q = ", Q, " done..."))
}


### Function to run regressions between classic metrics

god_BayReg_trait_met_ML2 <- function(resMetrics, Q, nMetrics = 5, pathSave,
                                     nChains, nIters, nCores, engine) {
  library(brms)
  library(cmdstanr)
  library(rstan)

  ### Inits
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  ### Store results
  fits <- list()
  R2_robust <- list()
  R2_regular <- list()

  for (i in 1:nMetrics) {
    print(names(resMetrics[i + 2]))
    print(names(resMetrics[i + 7]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 7)]
    headers <- names(data)

    formulas <- as.formula(paste(
      headers[3],
      " ~ ",
      paste(headers[4],
        paste0("+ (1|plotID)"),
        paste0("+ (1|Site)"),
        collapse = "+"
      )
    ))

    ### Run
    fit <- brms::brm(
      formula = formulas,
      data = data,
      iter = nIters,
      warmup = nIters / 5,
      cores = nCores,
      backend = engine
    )
    ## extract robust Bayes R2
    robust <- data.frame(bayes_R2(fit, robust = TRUE))
    robust$Trait <- headers[3]
    robust$Spec <- headers[4]
    R2_robust[[i]] <- robust

    ## extract Bayes R2
    reg <- data.frame(bayes_R2(fit))
    reg$Trait <- headers[3]
    reg$Spec <- headers[4]
    R2_regular[[i]] <- reg

    fits[[i]] <- fit
  }

  R2_robust <- do.call(rbind, R2_robust)
  R2_regular <- do.call(rbind, R2_regular)

  res <- list(fits = fits, R2_robust = R2_robust, R2 = R2_regular)

  save(res, file = pathSave)

  return(res)

  print(paste0("Bayesian models for Q = ", Q, " done..."))
}
