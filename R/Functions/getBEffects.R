
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