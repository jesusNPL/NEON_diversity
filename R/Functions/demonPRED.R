
### Function that correlates the observed and predited variables and
# calculates metrics of accuracy

demonPrediction <- function(fit, variable = "PD", covariable = "SD", dimension, Q,
                            groupLevel = FALSE, nCores = 12, nIters = 2000) {
  library(brms)
  library(cmdstanr)

  # Mape function
  # extracted from Jung 2022 PEERJ 10.7717/peerj.13872
  mape <- function(observed, predicted, type = "normal", denim = 1) {
    assertthat::assert_that(length(observed) == length(predicted))
    if (type == "normal") {
      return(
        mean(
          abs((observed - predicted) / observed),
          na.rm = TRUE
        ) * 100
      )
    } else {
      # Calculate symmetric map https://en.wikipedia.org/wiki/Symmetric_mean_absolute_percentage_error
      return(
        (100 / length(observed)) * sum(abs(predicted - observed) / ((abs(predicted) + abs(observed)) * denim),
          na.rm = TRUE
        )
      )
    }
  }

  rmse <- function(observed, predicted) {
    sqrt(mean((observed - predicted)^2))
  }

  ##### prepare data
  ## predicted values including the group-level effect of sites

  if (groupLevel == TRUE) {
    pp <- data.frame(predict(fit,
      robust = TRUE,
      probs = c(0.025, 0.05, 0.11, 0.25, 0.50, 0.75, 0.89, 0.95, 0.975)
    ))

    pp <- data.frame(fit$data, pp)

    headers <- names(pp)

    pp1 <- pp[, c(variable, "Estimate")]

    names(pp1) <- c("variable", "Estimate")
  } else {
    ## predicted values excluding the group-level effect of sites
    pp <- data.frame(predict(fit,
      re_formula = NA, #~ (1 | Site), 
      robust = TRUE,
      probs = c(0.025, 0.05, 0.11, 0.25, 0.50, 0.75, 0.89, 0.95, 0.975)
    ))

    pp <- data.frame(fit$data, pp)

    headers <- names(pp)

    pp1 <- pp[, c(variable, "Estimate")]

    names(pp1) <- c("variable", "Estimate")
  }

  # form <- as.formula(mvbind(paste(variable, "Estimate") ~ 1)))

  ## Correlation between observed and predicted
  corFit <- brm(
    data = pp1, family = student,
    formula = bf(mvbind("variable", "Estimate") ~ 1) +
      set_rescor(rescor = TRUE),
    iter = nIters, warmup = nIters / 5, chains = 4, cores = nCores,
    backend = "cmdstanr"
  )

  sumFit <- summary(corFit)
  cors <- round(sumFit$rescor_pars, 3)
  cors$variable <- variable
  cors$covariable <- covariable
  cors$predicted <- "Estimate"
  cors$dimension <- dimension
  cors$q <- Q
  ##### Metrics of accuracy and predictability all data #####

  ## Symetric Mean Absolute Percent Error
  smape <- mape(
    observed = pp[, variable],
    predicted = pp[, "Estimate"], type = "symetric"
  )

  ## RMSE - Root Mean Squared Error
  srmse <- rmse(observed = pp[, variable], predicted = pp[, "Estimate"])

  ## RMSLE - Root Mean Squared Log Error
  srmsle <- Metrics::rmsle(
    actual = pp[, variable],
    predicted = pp[, "Estimate"]
  )

  ## MAE - mean absolute error
  smae <- Metrics::mae(
    actual = pp[, variable],
    predicted = pp[, "Estimate"]
  )

  ## Bias
  sbias <- Metrics::bias(
    actual = pp[, variable],
    predicted = pp[, "Estimate"]
  )

  ### Add results to correlation data frame
  cors$smape <- smape
  cors$rmse <- srmse
  cors$rmlse <- srmsle
  cors$mae <- smae
  cors$bias <- sbias

  rownames(cors) <- NULL

  ##### Plot level evaluations #####
  PREDdev <- abs(pp[, variable] - pp[, "Estimate"]) / max(pp[, variable])
  PREDdiff <- (pp[, variable] - pp[, "Estimate"])
  PRED_change <- ((pp[, "Estimate"] - pp[, variable]) / pp[, variable])
  # negative PRED_change values suggest that community underprediction
  # positive PRED_change values suggest that community overprediction
  RRL <- log(pp[, "Estimate"] / pp[, variable])

  plotRES <- data.frame(pp,
    Deviation = PREDdev, Difference = PREDdiff,
    Change = PRED_change, LRR = RRL
  )

  plotRES$variable <- variable
  plotRES$covariable <- covariable
  plotRES$dimension <- dimension
  plotRES$q <- Q
  ### Store results
  results <- list(corTable = cors, siteResults = plotRES)

  return(results)
}

## Use
# load("Results/Regressions/phylo-spec/SAM/Scale/reg_PhyloSpec_DIS_SAM_scale_q0.RData")

# x <- demonPrediction(fit = res$fits[[1]], variable = "PD", covariable = "SD",
#                     dimension = "Phylo", nCores = 12, nIters = 2000)
