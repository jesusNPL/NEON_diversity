god_BayReg_alpha_taxo <- function(resMetrics, nMetrics, pathSave,
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
    print(names(resMetrics[i + 15]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 15)]
    headers <- names(data)

    formulas <- as.formula(paste(
      headers[3],
      " ~ ",
      paste(headers[4],
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
    robust$Taxo <- headers[3]
    robust$Spec <- headers[4]
    R2_robust[[i]] <- robust

    ## extract Bayes R2
    reg <- data.frame(bayes_R2(fit))
    reg$Taxo <- headers[3]
    reg$Spec <- headers[4]
    R2_regular[[i]] <- reg

    fits[[i]] <- fit
  }

  R2_robust <- do.call(rbind, R2_robust)
  R2_regular <- do.call(rbind, R2_regular)

  res <- list(fits = fits, R2_robust = R2_robust, R2 = R2_regular)

  save(res, file = pathSave)

  return(res)

  print(paste0("Bayesian models for alpha taxonomic dimension done..."))
}

############ ---------------- Alpha-Beta-Gamma ---------------- ##############
#q0
god_BayReg_ABG_taxo <- function(resMetrics, nMetrics, pathSave,
                                diversityLevel, qlevel,
                                nChains, nIters, nCores, engine) {
  library(brms)
  library(cmdstanr)
  library(rstan)
  library(tidyverse)
  
  ### Preparing data
  
  resMetrics <- resMetrics %>% 
    filter(Diversity == diversityLevel & q == qlevel)

  ### Inits
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  ### Store results
  fits <- list()
  R2_robust <- list()
  R2_regular <- list()

  for (i in 1:nMetrics) {
    print(names(resMetrics[i + 3]))
    print(names(resMetrics[i + 5]))

    ### Inits
    data <- resMetrics[, c(1:3, i + 3, i + 5)]
    headers <- names(data)

    formulas <- as.formula(paste(
      headers[4],
      " ~ ",
      paste(headers[5]),
        #paste0("+ (1|Site)"),
        collapse = "+"
      )
    )#)

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
    robust$Taxo <- headers[4]
    robust$Spec <- headers[5]
    R2_robust[[i]] <- robust

    ## extract Bayes R2
    reg <- data.frame(bayes_R2(fit))
    reg$Taxo <- headers[4]
    reg$Spec <- headers[5]
    R2_regular[[i]] <- reg

    fits[[i]] <- fit
  }

  R2_robust <- do.call(rbind, R2_robust)
  R2_regular <- do.call(rbind, R2_regular)

  res <- list(fits = fits, R2_robust = R2_robust, R2 = R2_regular)

  save(res, file = pathSave)

  return(res)

  print(paste0("Bayesian models for ",  diversityLevel, 
               " taxonomic dimension at level ", qlevel, " done..."))
}