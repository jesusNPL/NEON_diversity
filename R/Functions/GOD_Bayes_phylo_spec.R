god_BayReg_phylo <- function(resMetrics, Q, nMetrics, pathSave, scale = FALSE,
                             nChains, nIters, nCores, control, engine) {
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
    print(names(resMetrics[i + 12]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 12)]
    headers <- names(data)

    if (scale == FALSE) {
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
      control = control,
      backend = engine
    )
    ## extract robust Bayes R2
    robust <- data.frame(bayes_R2(fit, robust = TRUE))
    robust$Phylo <- headers[3]
    robust$Spec <- headers[4]
    R2_robust[[i]] <- robust

    ## extract Bayes R2
    reg <- data.frame(bayes_R2(fit))
    reg$Phylo <- headers[3]
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

############ ------------------- Multilevel 2 ------------------- ##############

god_BayReg_phylo_ML2 <- function(resMetrics, Q, nMetrics, pathSave,
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
    print(names(resMetrics[i + 12]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 12)]
    headers <- names(data)

    # formulas <- as.formula(paste(
    # headers[3],
    # " ~ ",
    # paste(headers[4],
    # paste0("+ (1|plotID)"),
    # paste0("+ (1|Site)"),
    # collapse = "+"
    # )
    # ))

    formulas <- as.formula(paste(headers[3],  " ~ ",
      paste0("scale(", headers[4], ")"),
      paste0("+ (1|plotID)"),
      paste0("+ (1|Site)"),
      collapse = "+"
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
    robust$Phylo <- headers[3]
    robust$Spec <- headers[4]
    R2_robust[[i]] <- robust

    ## extract Bayes R2
    reg <- data.frame(bayes_R2(fit))
    reg$Phylo <- headers[3]
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
