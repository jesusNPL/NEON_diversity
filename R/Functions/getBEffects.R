
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
  
  for(i in 1:length(fits)) { 
    
    print(estimates[i, 5])
    
    if (robust == TRUE) { 
      fix <- round(data.frame(brms::fixef(fits[[i]], 
                              robust = TRUE, 
                              probs = probs)), 3) 
    } else {
      fix <- round(data.frame(brms::fixef(fits[[i]], 
                              robust = FALSE, 
                              probs = probs)), 3)
      }
    
    fix$metric <- estimates[i, 5] 
    fix$params <- params
    fix$dimension <- dimension 
    fix$level <- level
    fix$Q <- Q
    fixLst[[i]] <- fix
    
  }
  fixdt <- do_call(rbind, fixLst) 
  rownames(fixdt) <- NULL
  return(fixdt)
}

#x <- getFixef(fits = res$fits, robust = TRUE, 
 #        estimates = res$R2_robust, 
  #       dimension = "phylogeny", level = "forest", Q = 0)

##### Function to get fixed effects from posterior distribution #####
getFixefPosterior <- function(fits, 
                              estimates, 
                              params = "^b", 
                              nDraws, 
                              dimension, 
                              Q, 
                              level) { 
  
  posteriorLST <- list()
  
  for(i in 1:length(fits)) {
    print(estimates[i, 5]) 
    posterior <- brms::posterior_samples(fits[[i]], pars = params) 
    keep <- sample(nrow(posterior), nDraws)
    posterior <- posterior[keep, ] 
    
    names(posterior) <- c("Incercept", "Slope")
    posterior$metric <- estimates[i, 5] 
    posterior$dimension <- dimension 
    posterior$level <- level
    posterior$Q <- Q
    
    posteriorLST[[i]] <- posterior
    
  }
  
  posteriorSamples <- do.call(rbind, posteriorLST)
  return(posteriorSamples)
}


#y <- getFixefPosterior(fits = res$fits, 
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

#getR2Estimates(fits = res, 
 #              robust = TRUE, 
  #             dimension = "phylo", 
   #            Q = 1, 
    #           level = "forest")

##### Function for Hypothesis testing #####

makeHypothesis <- function(fits, 
                           estimates, 
                           dimension, 
                           Q, 
                           level) { 
  
  hypLST <- list()
  
  for(i in 1:length(fits)) {
    
    fit <- fits[[i]]
    
    fixfit <- data.frame(fixef(fit))
    param <- rownames(fixfit)[2]
    slope <- fixfit[2, 1]
    
    if (slope > 0) {
      hyp <- paste0(param, " > 0") 
      h <- hypothesis(fit, hyp, class = "b")$hypothesis
    } else {
      hyp <- paste0(param, " < 0") 
      h <- hypothesis(fit, hyp, class = "b")$hypothesis
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

#x <- makeHypothesis(fits = res$fits, 
 #                   estimates = res$R2_robust, 
  #                  dimension = "phylogeny", 
   #                 Q = "q2", 
    #                level = "all")
