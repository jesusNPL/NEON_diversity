###### Bayesian moments #####
god_BayReg_moments <- function(momentsRES, nChains, nIters, nCores, engine, nMetrics,
                               pathSaveTaxo, pathSaveTrait, pathSavePhylo) {
  library(brms)
  library(rstan)
  library(cmdstanr)

  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  # Set prior for Bayesian robust correlations
  Corr_prior <- c(
    prior(gamma(2, 0.1), class = nu),
    prior(lkj(1), class = rescor)
  )

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

  for (i in 1:nMetrics) {
    print("Taxonomy ~ Spectra")
    print(names(taxo[i + 2]))
    print(names(taxo[i + 8]))

    dataTx <- taxo[, c(1:2, i + 2, i + 8)]
    headersTx <- names(dataTx)
    formTx <- as.formula(paste(
      headersTx[3], " ~ ",
      paste(headersTx[4],
        paste0("+ (1|Site)"),
        collapse = "+"
      )
    ))

    fitTx <- brm(formTx,
      data = dataTx,
      # family = student(),
      iter = nIters,
      warmup = nIters / 5,
      chains = nChains,
      cores = nCores,
      backend = engine
    )

    robustTaxo <- data.frame(bayes_R2(fitTx, robust = TRUE))
    robustTaxo$Taxo <- headersTx[3]
    robustTaxo$Spec <- headersTx[4]
    R2_robust_taxo[[i]] <- robustTaxo

    regularTaxo <- data.frame(bayes_R2(fitTx))
    regularTaxo$Taxo <- headersTx[3]
    regularTaxo$Spec <- headersTx[4]
    R2_regular_taxo[[i]] <- regularTaxo

    fits_taxo[[i]] <- fitTx

    print("Trait ~ Spectra")
    print(names(trait[i + 2]))
    print(names(trait[i + 8]))

    dataTrt <- trait[, c(1:2, i + 2, i + 8)]
    headersTrt <- names(dataTrt)
    formTrt <- as.formula(paste(
      headersTrt[3], " ~ ",
      paste(headersTrt[4],
        paste0("+ (1|Site)"),
        collapse = "+"
      )
    ))

    fitTrt <- brm(formTrt,
      data = dataTrt,
      # family = student(),
      iter = nIters,
      warmup = nIters / 5,
      chains = nChains,
      cores = nCores,
      backend = engine
    )

    robustTrait <- data.frame(bayes_R2(fitTrt, robust = TRUE))
    robustTrait$Trait <- headersTrt[3]
    robustTrait$Spec <- headersTrt[4]
    R2_robust_trait[[i]] <- robustTrait

    regularTrait <- data.frame(bayes_R2(fitTrt))
    regularTrait$Taxo <- headersTrt[3]
    regularTrait$Spec <- headersTrt[4]
    R2_regular_trait[[i]] <- regularTrait

    fits_trait[[i]] <- fitTrt

    print("Phylogeny ~ Spectra")
    print(names(phylo[i + 2]))
    print(names(phylo[i + 8]))

    dataPhy <- phylo[, c(1:2, i + 2, i + 8)]
    headersPhy <- names(dataPhy)
    formPhy <- as.formula(paste(
      headersPhy[3], " ~ ",
      paste(headersPhy[4],
        paste0("+ (1|Site)"),
        collapse = "+"
      )
    ))

    fitPhy <- brm(formPhy,
      data = dataPhy,
      # family = student(),
      iter = nIters,
      warmup = nIters / 5,
      chains = nChains,
      cores = nCores,
      backend = engine
    )

    robustPhylo <- data.frame(bayes_R2(fitPhy, robust = TRUE))
    robustPhylo$Phylo <- headersPhy[3]
    robustPhylo$Spec <- headersPhy[4]
    R2_robust_phylo[[i]] <- robustPhylo

    regularPhylo <- data.frame(bayes_R2(fitPhy))
    regularPhylo$Phylo <- headersPhy[3]
    regularPhylo$Spec <- headersPhy[4]
    R2_regular_phylo[[i]] <- regularPhylo

    fits_phylo[[i]] <- fitPhy
  }
  # Taxonomy
  R2_robust_taxo <- do.call(rbind, R2_robust_taxo)
  R2_regular_taxo <- do.call(rbind, R2_regular_taxo)

  resultsTaxo <- list(
    CorrRobust = R2_robust_taxo,
    Corr = R2_regular_taxo,
    fits = fits_taxo
  )
  save(resultsTaxo, file = pathSaveTaxo)
  # Trait
  R2_robust_trait <- do.call(rbind, R2_robust_trait)
  R2_regular_trait <- do.call(rbind, R2_regular_trait)

  resultsTrait <- list(
    CorrRobust = R2_robust_trait,
    Corr = R2_regular_trait,
    fits = fits_trait
  )
  save(resultsTrait, file = pathSaveTrait)
  # Phylogeny
  R2_robust_phylo <- do.call(rbind, R2_robust_phylo)
  R2_regular_phylo <- do.call(rbind, R2_regular_phylo)

  resultsPhylo <- list(
    CorrRobust = R2_robust_phylo,
    Corr = R2_regular_phylo,
    fits = fits_phylo
  )
  save(resultsPhylo, file = pathSavePhylo)
  # return(results)

  print("Bayesian models for multiple moments... done!!!")
}

########## ------------- Multilevel 2 ----------------- ###############

###### Bayesian moments #####
god_BayReg_moments_ML2 <- function(momentsRES, nChains, nIters, nCores, engine, nMetrics,
                                   pathSaveTaxo, pathSaveTrait, pathSavePhylo) {
  library(brms)
  library(rstan)
  library(cmdstanr)

  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  # Set prior for Bayesian robust correlations
  Corr_prior <- c(
    prior(gamma(2, 0.1), class = nu),
    prior(lkj(1), class = rescor)
  )

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

  for (i in 1:nMetrics) {
    print("Taxonomy ~ Spectra")
    print(names(taxo[i + 2]))
    print(names(taxo[i + 8]))

    dataTx <- taxo[, c(1:2, i + 2, i + 8)]
    headersTx <- names(dataTx)
    formTx <- as.formula(paste(
      headersTx[3], " ~ ",
      paste(headersTx[4],
        paste0("+ (1|plotID)"),
        paste0("+ (1|Site)"),
        collapse = "+"
      )
    ))

    fitTx <- brm(formTx,
      data = dataTx,
      # family = student(),
      iter = nIters,
      warmup = nIters / 5,
      chains = nChains,
      cores = nCores,
      backend = engine
    )

    robustTaxo <- data.frame(bayes_R2(fitTx, robust = TRUE))
    robustTaxo$Taxo <- headersTx[3]
    robustTaxo$Spec <- headersTx[4]
    R2_robust_taxo[[i]] <- robustTaxo

    regularTaxo <- data.frame(bayes_R2(fitTx))
    regularTaxo$Taxo <- headersTx[3]
    regularTaxo$Spec <- headersTx[4]
    R2_regular_taxo[[i]] <- regularTaxo

    fits_taxo[[i]] <- fitTx

    print("Trait ~ Spectra")
    print(names(trait[i + 2]))
    print(names(trait[i + 8]))

    dataTrt <- trait[, c(1:2, i + 2, i + 8)]
    headersTrt <- names(dataTrt)
    formTrt <- as.formula(paste(
      headersTrt[3], " ~ ",
      paste(headersTrt[4],
        paste0("+ (1|plotID)"),
        paste0("+ (1|Site)"),
        collapse = "+"
      )
    ))

    fitTrt <- brm(formTrt,
      data = dataTrt,
      # family = student(),
      iter = nIters,
      warmup = nIters / 5,
      chains = nChains,
      cores = nCores,
      backend = engine
    )

    robustTrait <- data.frame(bayes_R2(fitTrt, robust = TRUE))
    robustTrait$Trait <- headersTrt[3]
    robustTrait$Spec <- headersTrt[4]
    R2_robust_trait[[i]] <- robustTrait

    regularTrait <- data.frame(bayes_R2(fitTrt))
    regularTrait$Taxo <- headersTrt[3]
    regularTrait$Spec <- headersTrt[4]
    R2_regular_trait[[i]] <- regularTrait

    fits_trait[[i]] <- fitTrt

    print("Phylogeny ~ Spectra")
    print(names(phylo[i + 2]))
    print(names(phylo[i + 8]))

    dataPhy <- phylo[, c(1:2, i + 2, i + 8)]
    headersPhy <- names(dataPhy)
    formPhy <- as.formula(paste(
      headersPhy[3], " ~ ",
      paste(headersPhy[4],
        paste0("+ (1|plotID)"),
        paste0("+ (1|Site)"),
        collapse = "+"
      )
    ))

    fitPhy <- brm(formPhy,
      data = dataPhy,
      # family = student(),
      iter = nIters,
      warmup = nIters / 5,
      chains = nChains,
      cores = nCores,
      backend = engine
    )

    robustPhylo <- data.frame(bayes_R2(fitPhy, robust = TRUE))
    robustPhylo$Phylo <- headersPhy[3]
    robustPhylo$Spec <- headersPhy[4]
    R2_robust_phylo[[i]] <- robustPhylo

    regularPhylo <- data.frame(bayes_R2(fitPhy))
    regularPhylo$Phylo <- headersPhy[3]
    regularPhylo$Spec <- headersPhy[4]
    R2_regular_phylo[[i]] <- regularPhylo

    fits_phylo[[i]] <- fitPhy
  }
  # Taxonomy
  R2_robust_taxo <- do.call(rbind, R2_robust_taxo)
  R2_regular_taxo <- do.call(rbind, R2_regular_taxo)

  resultsTaxo <- list(
    CorrRobust = R2_robust_taxo,
    Corr = R2_regular_taxo,
    fits = fits_taxo
  )
  save(resultsTaxo, file = pathSaveTaxo)
  # Trait
  R2_robust_trait <- do.call(rbind, R2_robust_trait)
  R2_regular_trait <- do.call(rbind, R2_regular_trait)

  resultsTrait <- list(
    CorrRobust = R2_robust_trait,
    Corr = R2_regular_trait,
    fits = fits_trait
  )
  save(resultsTrait, file = pathSaveTrait)
  # Phylogeny
  R2_robust_phylo <- do.call(rbind, R2_robust_phylo)
  R2_regular_phylo <- do.call(rbind, R2_regular_phylo)

  resultsPhylo <- list(
    CorrRobust = R2_robust_phylo,
    Corr = R2_regular_phylo,
    fits = fits_phylo
  )
  save(resultsPhylo, file = pathSavePhylo)
  # return(results)

  print("Bayesian models for multiple moments... done!!!")
}
