##### Function to run regressions phylogenetic dimension #####
god_BayReg_local_phylo <- function(resMetrics, Q, nMetrics, scale = FALSE,
                                   site, dimension, nChains, nIters, nCores,
                                   control, engine, alpha = 0.05) {
  library(brms)
  library(cmdstanr)
  library(rstan)

  ### Inits
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  ### Store results
  R2_robust <- list()

  fixef_robust <- list()

  hypLST <- list()

  for (i in 1:nMetrics) {
    print(names(resMetrics[i + 2]))
    print(names(resMetrics[i + 9]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 9)]
    headers <- names(data)

    if (scale == FALSE) {
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste0("t2(", headers[4], ")")
      ))
    } else {
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste0("scale(", headers[4], ")")
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
    robust$Site <- site
    R2_robust[[i]] <- robust

    ## extract Fixef
    probs <- c(0.025, 0.05, 0.11, 0.25, 0.75, 0.89, 0.95, 0.975)

    fix <- round(data.frame(brms::fixef(fit,
      robust = TRUE,
      probs = probs
    )), 3)
    fix$Phylo <- headers[3]
    fix$Spec <- headers[4]
    fix$formula <- paste(headers[3], " ~ ", headers[4])
    fix$Site <- site
    fixef_robust[[i]] <- fix

    ## Make hypothesis
    param <- rownames(fix)[2]
    slope <- fix[2, 1]

    if (slope > 0) {
      hyp <- paste0(param, " > 0")
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    } else {
      hyp <- paste0(param, " < 0")
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    }

    h$metric <- paste(headers[3], " ~ ", headers[4])
    h$dimension <- dimension
    h$Q <- Q
    h$Site <- site
    h$alpha <- alpha
    
    hypLST[[i]] <- h
  }

  R2_robust <- do.call(rbind, R2_robust)
  fixef_robust <- do.call(rbind, fixef_robust)
  hypTBL <- do.call(rbind, hypLST)

  res <- list(
    R2_robust = R2_robust,
    fixef_robust = fixef_robust,
    hypTBL = hypTBL
  )

  return(res)

  print(paste0("Bayesian models for site", site, " at Q = ", Q, " done..."))
}

##### Function to run regressions between simple distance trait metrics #####
god_BayReg_trait_local_dis <- function(resMetrics, Q, nMetrics,
                                       scale = FALSE, site, dimension,
                                       nChains, nIters, nCores, engine, alpha = 0.05) {
  library(brms)
  library(cmdstanr)
  library(rstan)

  ### Inits
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  ### Store results
  R2_robust <- list()
  fixef_robust <- list()
  hypLST <- list()

  for (i in 1:nMetrics) {
    print(names(resMetrics[i + 2]))
    print(names(resMetrics[i + 8]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 8)]
    headers <- names(data)

    if (scale == FALSE) {
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste0("t2(", headers[4], ")")
      ))
    } else {
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste0("scale(", headers[4], ")")
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
    robust$Phylo <- headers[3]
    robust$Spec <- headers[4]
    robust$Site <- site
    R2_robust[[i]] <- robust

    ## extract Fixef
    probs <- c(0.025, 0.05, 0.11, 0.25, 0.75, 0.89, 0.95, 0.975)

    fix <- round(data.frame(brms::fixef(fit,
      robust = TRUE,
      probs = probs
    )), 3)

    fix$Phylo <- headers[3]
    fix$Spec <- headers[4]
    fix$formula <- paste(headers[3], " ~ ", headers[4])
    fix$Site <- site
    fixef_robust[[i]] <- fix

    ## Make hypothesis
    param <- rownames(fix)[2]
    slope <- fix[2, 1]

    if (slope > 0) {
      hyp <- paste0(param, " > 0")
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    } else {
      hyp <- paste0(param, " < 0")
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    }

    h$metric <- paste(headers[3], " ~ ", headers[4])
    h$dimension <- dimension
    h$Q <- Q
    h$Site <- site 
    h$alpha <- alpha

    hypLST[[i]] <- h
  }

  R2_robust <- do.call(rbind, R2_robust)
  fixef_robust <- do.call(rbind, fixef_robust)
  hypTBL <- do.call(rbind, hypLST)

  res <- list(
    R2_robust = R2_robust,
    fixef_robust = fixef_robust,
    hypTBL = hypTBL
  )

  return(res)

  print(paste0("Bayesian models for site", site, " at Q = ", Q, " done..."))
}

##### Function to run regressions between classic trait metrics #####
god_BayReg_trait_local_met <- function(resMetrics, Q, nMetrics = 5,
                                       scale = FALSE, site, dimension,
                                       nChains, nIters, nCores, engine, alpha = 0.05) {
  library(brms)
  library(cmdstanr)
  library(rstan)

  ### Inits
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  ### Store results
  R2_robust <- list()
  fixef_robust <- list()
  hypLST <- list()

  for (i in 1:nMetrics) {
    print(names(resMetrics[i + 2]))
    print(names(resMetrics[i + 7]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 7)]
    headers <- names(data)

    if (scale == FALSE) {
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste0("t2(", headers[4], ")")
      ))
    } else {
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste0("scale(", headers[4], ")")
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
    robust$Phylo <- headers[3]
    robust$Spec <- headers[4]
    robust$Site <- site
    R2_robust[[i]] <- robust

    ## extract Fixef
    probs <- c(0.025, 0.05, 0.11, 0.25, 0.75, 0.89, 0.95, 0.975)

    fix <- round(data.frame(brms::fixef(fit,
      robust = TRUE,
      probs = probs
    )), 3)

    fix$Phylo <- headers[3]
    fix$Spec <- headers[4]
    fix$formula <- paste(headers[3], " ~ ", headers[4])
    fix$Site <- site
    fixef_robust[[i]] <- fix

    ## Make hypothesis
    param <- rownames(fix)[2]
    slope <- fix[2, 1]

    if (slope > 0) {
      hyp <- paste0(param, " > 0")
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    } else {
      hyp <- paste0(param, " < 0")
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    }

    h$metric <- paste(headers[3], " ~ ", headers[4])
    h$dimension <- dimension
    h$Q <- Q
    h$Site <- site
    h$alpha <- alpha

    hypLST[[i]] <- h
  }

  R2_robust <- do.call(rbind, R2_robust)
  fixef_robust <- do.call(rbind, fixef_robust)
  hypTBL <- do.call(rbind, hypLST)

  res <- list(
    R2_robust = R2_robust,
    fixef_robust = fixef_robust,
    hypTBL = hypTBL
  )

  return(res)

  print(paste0("Bayesian models for site", site, " at Q = ", Q, " done..."))
}

##### Function to run regressions taxonomic dimension ######
god_BayReg_alpha_local_taxo <- function(resMetrics, nMetrics, scale = FALSE,
                                        site, dimension,
                                        nChains, nIters, nCores, engine, alpha = 0.05) {
  library(brms)
  library(cmdstanr)
  library(rstan)

  ### Inits
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  ### Store results
  R2_robust <- list()
  fixef_robust <- list()
  hypLST <- list()

  for (i in 1:nMetrics) {
    print(names(resMetrics[i + 2]))
    print(names(resMetrics[i + 15]))

    ### Inits
    data <- resMetrics[, c(1:2, i + 2, i + 15)]
    headers <- names(data)

    if (scale == FALSE) {
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste0("t2(", headers[4], ")")
      ))
    } else {
      formulas <- as.formula(paste(
        headers[3], " ~ ",
        paste0("scale(", headers[4], ")")
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
    robust$Phylo <- headers[3]
    robust$Spec <- headers[4]
    robust$Site <- site
    R2_robust[[i]] <- robust

    ## extract Fixef
    probs <- c(0.025, 0.05, 0.11, 0.25, 0.75, 0.89, 0.95, 0.975)

    fix <- round(data.frame(brms::fixef(fit,
      robust = TRUE,
      probs = probs
    )), 3)

    fix$Phylo <- headers[3]
    fix$Spec <- headers[4]
    fix$formula <- paste(headers[3], " ~ ", headers[4])
    fix$Site <- site
    fixef_robust[[i]] <- fix

    ## Make hypothesis
    param <- rownames(fix)[2]
    slope <- fix[2, 1]

    if (slope > 0) {
      hyp <- paste0(param, " > 0")
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    } else {
      hyp <- paste0(param, " < 0")
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    }

    h$metric <- paste(headers[3], " ~ ", headers[4])
    h$dimension <- dimension
    # h$Q <- Q
    h$Site <- site
    h$alpha <- alpha
    
    hypLST[[i]] <- h
  }

  R2_robust <- do.call(rbind, R2_robust)
  fixef_robust <- do.call(rbind, fixef_robust)
  hypTBL <- do.call(rbind, hypLST)

  res <- list(
    R2_robust = R2_robust,
    fixef_robust = fixef_robust,
    hypTBL = hypTBL
  )

  return(res)

  print(paste0("Bayesian models for site", site, " at Q = ", Q, " done..."))
}
