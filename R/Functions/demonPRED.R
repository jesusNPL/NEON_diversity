##### Function  demonPrediction() #####
# This function correlates the observed and predicted variables and calculates metrics of accuracy
# Fitted models are used to obtain predicted metric values

demonPrediction <- function(fit, 
                            variable = "PD", 
                            covariable = "SD", 
                            dimension, 
                            Q,
                            groupLevel = TRUE, 
                            nCores = 12, 
                            nIters = 2000, 
                            burnin) { 
  
  ### Main library
  library(brms)

  ### Mape function
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

  ### RSME function
  rmse <- function(observed, predicted) {
    sqrt(mean((observed - predicted)^2))
  }

  ##### prepare data #####
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
      re_formula = NA, 
      robust = TRUE,
      probs = c(0.025, 0.05, 0.11, 0.25, 0.50, 0.75, 0.89, 0.95, 0.975)
    ))

    pp <- data.frame(fit$data, pp)

    headers <- names(pp)

    pp1 <- pp[, c(variable, "Estimate")]

    names(pp1) <- c("variable", "Estimate") 
    
  }

  ### Correlation between observed and predicted metric values 
  corFit <- brm(
    data = pp1, 
    family = student,
    formula = bf(mvbind("variable", "Estimate") ~ 1) +
      set_rescor(rescor = TRUE),
    iter = nIters, 
    warmup = burning, 
    chains = 4, 
    cores = nCores,
    backend = "cmdstanr"
  )

  ### Extract correlations
  sumFit <- summary(corFit)
  cors <- round(sumFit$rescor_pars, 3)
  cors$variable <- variable
  cors$covariable <- covariable
  cors$predicted <- "Estimate"
  cors$dimension <- dimension
  cors$q <- Q
  
  ##### Metrics of accuracy and predictability all data #####
  ### Symetric Mean Absolute Percent Error
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
  
  ### Prediction deviation
  PREDdev <- abs(pp[, variable] - pp[, "Estimate"]) / max(pp[, variable])
  
  ### Prediction difference
  PREDdiff <- (pp[, variable] - pp[, "Estimate"])
  
  ### Prediction rate of change
  PRED_change <- ((pp[, "Estimate"] - pp[, variable]) / pp[, variable])
  # negative PRED_change values suggest that community underprediction
  # positive PRED_change values suggest that community overprediction
  
  ### Log ratio
  RRL <- log(pp[, "Estimate"] / pp[, variable])

  ### Store results at plot level
  plotRES <- data.frame(pp,
    Deviation = PREDdev, 
    Difference = PREDdiff,
    Change = PRED_change, 
    LRR = RRL
  )

  ## Add descriptions
  plotRES$variable <- variable
  plotRES$covariable <- covariable
  plotRES$dimension <- dimension
  plotRES$q <- Q 
  
  ### Store results
  results <- list(corTable = cors, siteResults = plotRES)

  return(results) 
  
} # end function 

## Usage
#load("Results/Regressions/Reanalyses/phylo_DIS/BLM_phylo_DIS_alpha_q0_scale.RData")

#x <- demonPrediction(fit = res$fits[[1]], 
 #                    variable = "PD", 
  #                   covariable = "SD",
   #                  dimension = "Phylogeny", 
                    # Q = Q, 
                    # groupLevel = TRUE, 
    #                 nCores = 12, 
     #                nIters = nIters, 
                    # burnin = nIters/5)


##### Function godPrediction() #####

# This function uses data partition to fit a model and then evaluates the model predictions

godPrediction <- function(trainData, # data to train the model 
                          testData, # data to test the model
                          variable = "PD", 
                          covariable = "SD", 
                          dimension, 
                          Q,
                          nCores, 
                          nIters, 
                          burnin, 
                          nChains, 
                          control, 
                          engine) { 
  
  ### Main library
  library(brms)
  
  ### Mape function
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
  
  ### RSME function
  rmse <- function(observed, predicted) {
    sqrt(mean((observed - predicted)^2))
  }
  
  ##### Model fitting using train data #####
  print(paste0("Model fitting using train data for ", variable, " ~ ", covariable))
  
  trainData <- trainData[, c("siteID", variable, covariable)] 
  names(trainData) <- c("siteID", "variable", "covariable")
  
  fit <- brm(
    formula = variable ~ covariable + (1|siteID), 
    data = trainData,
    iter = nIters,
    warmup = burnin,
    cores = nCores, 
    chains = nChains, 
    control = control,
    backend = engine
    ) 
  
  ##### Model predictions using test data #####  
  print(paste0("Model prediction using test data for ", variable, " ~ ", covariable))
  
  testData <- testData[, c("siteID", variable, covariable)] 
  names(testData) <- c("siteID", "variable", "covariable")
  
  pred <- predict(object = fit, 
                  newdata = testData, 
                  robust = TRUE, 
                  probs = c(0.025, 0.05, 0.11, 0.25, 0.50, 0.75, 0.89, 0.95, 0.975)) 
  
  pred <- data.frame(pred)
  
  ##### prepare data for correlations #####
  pp <- data.frame(testData, pred)
  ## predicted values including the group-level effect of sites 
  headers <- names(pp)
  
  pp1 <- pp[, c("variable", "Estimate")]
  
  names(pp1) <- c("variable", "Estimate") 
  
  ### Correlation between observed and predicted metric values 
  print("Correlation between predicted and observed metrics values")
  
  corFit <- brm(
    data = pp1, 
    family = student,
    formula = bf(mvbind("variable", "Estimate") ~ 1) +
      set_rescor(rescor = TRUE),
    iter = nIters, 
    warmup = burnin, 
    chains = nChains, 
    cores = nCores,
    backend = "cmdstanr"
  )
  
  ### Extract correlations
  sumFit <- summary(corFit)
  cors <- round(sumFit$rescor_pars, 3)
  cors$variable <- variable
  cors$covariable <- covariable
  cors$predicted <- "Estimate"
  cors$dimension <- dimension
  cors$Q <- Q
  
  ##### Metrics of accuracy and predictability all data #####
  
  ### Symetric Mean Absolute Percent Error
  smape <- mape(
    observed = pp[, "variable"],
    predicted = pp[, "Estimate"], 
    type = "symetric"
  )
  
  ## RMSE - Root Mean Squared Error
  srmse <- rmse(observed = pp[, "variable"], 
                predicted = pp[, "Estimate"])
  
  ## RMSLE - Root Mean Squared Log Error
  srmsle <- Metrics::rmsle(
    actual = pp[, "variable"],
    predicted = pp[, "Estimate"]
  )
  
  ## MAE - mean absolute error
  smae <- Metrics::mae(
    actual = pp[, "variable"],
    predicted = pp[, "Estimate"]
  )
  
  ## Bias
  sbias <- Metrics::bias(
    actual = pp[, "variable"],
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
  
  ### Prediction deviation
  PREDdev <- abs(pp[, "variable"] - pp[, "Estimate"]) / max(pp[, "variable"])
  
  ### Prediction difference
  PREDdiff <- (pp[, "variable"] - pp[, "Estimate"])
  
  ### Prediction rate of change
  PRED_change <- ((pp[, "Estimate"] - pp[, "variable"]) / pp[, "variable"])
  # negative PRED_change values suggest that community underprediction
  # positive PRED_change values suggest that community overprediction
  
  ### Log ratio
  RRL <- log(pp[, "Estimate"] / pp[, "variable"])
  
  ### Store results at plot level
  plotRES <- data.frame(pp,
                        Deviation = PREDdev, 
                        Difference = PREDdiff,
                        Change = PRED_change, 
                        LRR = RRL
  )
  
  ## Add descriptions
  plotRES$Ground_metric <- variable
  plotRES$Spectral_metric <- covariable
  plotRES$dimension <- dimension
  plotRES$Q <- Q 
  
  ### Reorder results at plot level
  plotRES <- plotRES %>% 
    select(siteID, Ground_metric, Spectral_metric, dimension, Q, everything())
  
  ### Store results
  results <- list(corTable = cors, siteResults = plotRES)
  
  return(results) 
  
} # end function
