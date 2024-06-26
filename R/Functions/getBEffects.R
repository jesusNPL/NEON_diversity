
##### Function to extract fixed effects and credible intervals #####
getFixef <- function(fits,
                     probs = c(0.025, 0.05, 0.11, 0.25, 0.75, 0.89, 0.95, 0.975),
                     robust = TRUE,
                     estimates,
                     params = c("Intercept", "Slope"),
                     dimension,
                     Q,
                     level) { 
  
  fixLst <- list()

  for (i in 1:length(fits)) { 
    
    print(estimates[i, 5])

    if (robust == TRUE) { 
      
      fix <- round(data.frame(brms::fixef(fits[[i]],
        robust = TRUE,
        probs = probs
      )), 3) 
      
    } else { 
      
      fix <- round(data.frame(brms::fixef(fits[[i]],
        robust = FALSE,
        probs = probs
      )), 3) 
      
    }

    fix$metric <- estimates[i, 5]
    fix$params <- params
    fix$dimension <- dimension
    fix$level <- level
    fix$Q <- Q
    fixLst[[i]] <- fix 
    
  } 
  
  fixdt <- do.call(rbind, fixLst)
  rownames(fixdt) <- NULL 
  
  fixdt <- fixdt %>% 
    dplyr::select(dimension, metric, params, level, Q, everything())
  
  return(fixdt) 
  
}

# x <- getFixef(fits = res$fits, robust = TRUE,
#        estimates = res$R2_robust,
#       dimension = "phylogeny", level = "forest", Q = 0)

##### Function to extract fixed effects and credible intervals #####
getFixef_simple <- function(fit,
                     probs = c(0.025, 0.05, 0.11, 0.25, 0.75, 0.89, 0.95, 0.975),
                     robust = TRUE,
                     metric,
                     estimates,
                     params = c("Intercept", "Slope"),
                     dimension,
                     Q,
                     level) {
  
    
    if (robust == TRUE) { 
      
      fix <- round(data.frame(brms::fixef(fit,
                                          robust = TRUE,
                                          probs = probs
      )), 3) 
      
    } else { 
      
      fix <- round(data.frame(brms::fixef(fit,
                                          robust = FALSE,
                                          probs = probs
      )), 3) 
      
    }
    
    fix$metric <- metric
    fix$params <- params
    fix$dimension <- dimension
    fix$Quantile <- level
    fix$Q <- Q

  fixdt <- fix
  rownames(fixdt) <- NULL
  return(fixdt) 
  
}

##### Function to get fixed effects from posterior distribution #####
getFixefPosterior <- function(fits,
                              estimates,
                              params = "^b",
                              nDraws,
                              dimension,
                              Q,
                              level) { 
  
  posteriorLST <- list()

  for (i in 1:length(fits)) { 
    
    print(estimates[i, 5]) 
    
    #posterior <- brms::posterior_samples(fits[[i]], pars = params) 
    ## Samples from posterior distribution
    posterior <- as.data.frame(x = fits[[i]], # fit object
                               regex = TRUE) # transform the parameters to a data frame
    
    ## Select samples to keep
    keep <- sample(nrow(posterior), nDraws)
    ## Keep samples at random
    posterior <- posterior[keep, ]
    ## Select columns for our purpose
    posterior <- posterior[, ]
    names(posterior) <- c("Incercept", "Slope", "Sigma") 
    
    posterior$metric <- estimates[i, 5] 
    posterior$dimension <- dimension 
    posterior$level <- level 
    posterior$Q <- Q 

    posteriorLST[[i]] <- posterior 
    
  }
  
  # Combine samples
  posteriorSamples <- do.call(rbind, posteriorLST) 
  
  return(posteriorSamples) 
  
}


# y <- getFixefPosterior(fits = res$fits,
#                      estimates = res$R2_robust,
#                     nDraws = 10,
#                    dimension = "phylo",
#                   Q = 0,
#                  level = "forest")

##### Function to get R2 from models #####
getR2Estimates <- function(fits,
                           robust = TRUE,
                           dimension,
                           Q,
                           level) { 
  
  if (robust == TRUE) { 
    
    estimates <- res$R2_robust 
    
  } else { 
    
    estimates <- res$R2 
    
  }

  names(estimates) <- c("Estimate", "Est.Error", "Q2.5", "Q97.5", "metric", "metric_spec")
  estimates$dimension <- dimension
  estimates$level <- level
  estimates$Q <- Q

  rownames(estimates) <- NULL
  return(estimates) 
  
}

# getR2Estimates(fits = res,
#              robust = TRUE,
#             dimension = "phylo",
#            Q = 1,
#           level = "forest")

##### Function for Hypothesis testing #####

makeHypothesis <- function(fits,
                           estimates,
                           dimension,
                           Q,
                           level, 
                           alpha = 0.05) { 
  
  hypLST <- list()

  for (i in 1:length(fits)) { 
    
    fit <- fits[[i]]

    fixfit <- data.frame(fixef(fit))
    param <- rownames(fixfit)[2]
    slope <- fixfit[2, 1]

    if (slope > 0) { 
      
      hyp <- paste0(param, " > 0") 
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis 
      
    } else { 
      
      hyp <- paste0(param, " < 0") 
      h <- hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis 
      
    }

    h$metric <- estimates[i, 5]
    h$dimension <- dimension
    h$level <- level
    h$Q <- Q

    hypLST[[i]] <- h 
    
  }

  hypTBL <- do.call(rbind, hypLST)
  return(hypTBL) 
  
}

# x <- makeHypothesis(fits = res$fits,
#                   estimates = res$R2_robust,
#                  dimension = "phylogeny",
#                 Q = "q2",
#                level = "all")

##### Function to make tables from fixed effects site level #####

makeTableSite <- function(resSite, Q, dimension) {
  params <- c("Intercept", "Slope")

  r2 <- list()
  slopes <- list()
  hypothesis <- list()

  for (j in 1:length(resSite)) {
    esti <- resSite[[j]]$R2_robust
    esti$Q <- Q
    esti$dimension <- dimension
    r2[[j]] <- esti

    betas <- resSite[[j]]$fixef_robust
    betas$params <- rep(params, nrow(resSite[[j]]$fixef_robust) / 2)
    betas$Q <- Q
    betas$dimension <- dimension
    slopes[[j]] <- betas

    hypothesis[[j]] <- resSite[[j]]$hypTBL
  }

  r2 <- do.call(rbind, r2)
  rownames(r2) <- NULL
  r2 <- r2 %>%
    select(Site, Q, Phylo, Spec, everything())

  slopes <- do.call(rbind, slopes)
  rownames(slopes) <- NULL
  slopes <- slopes %>%
    select(Site, Q, formula, params, everything())

  hypothesis <- do.call(rbind, hypothesis)
  rownames(hypothesis) <- NULL
  hypothesis <- hypothesis %>% 
    mutate(Q = Q) %>% 
    select(Site, Q, metric, everything())

  res <- list(estimate = r2, betas = slopes, hypothesis = hypothesis)
  return(res)
}

# x <- makeTableSite(resSite = res_phylo_q0, Q = "q0", dimension = "phylogeny")

##### Function to extract fixed effects and credible intervals #####
getFixefQuantREG <- function(fits,
                             probs = c(0.025, 0.05, 0.11, 0.25, 0.75, 0.89, 0.95, 0.975),
                             robust = TRUE,
                             metric,
                             params = c("Intercept", "Slope"),
                             dimension,
                             Q,
                             QuantLevel) { 
  
  fixLst <- list()
  
  for (i in 1:length(fits)) { 
    
    print(paste(metric, QuantLevel[i]))
    
    if (robust == TRUE) { 
      
      fix <- round(data.frame(brms::fixef(fits[[i]],
                                          robust = TRUE,
                                          probs = probs
      )), 3) 
      
    } else { 
      
      fix <- round(data.frame(brms::fixef(fits[[i]],
                                          robust = FALSE,
                                          probs = probs
      )), 3) 
      
    }
    
    fix$metric <- metric
    fix$params <- params
    fix$dimension <- dimension
    fix$Quantile <- QuantLevel[i]
    fix$Q <- Q
    fixLst[[i]] <- fix 
    
  } 
  
  fixdt <- do.call(rbind, fixLst)
  rownames(fixdt) <- NULL
  return(fixdt) 
  
}

