library(brms)
library("bayesplot")
theme_set(bayesplot::theme_default(base_family = "sans"))

# to extract group-level standard deviations use "^sd_"
# to extract posterior of population-level effect use "^b"

## Auxiliary function that allows extract parameters from Bayesian models
extractSamples <- function(model, parameters = "^b", nDraws = 100) { 
  
  #posterior <- as.matrix(model, pars = parameters)
  posterior <- posterior_samples(model, pars = parameters)
  post_means <- colMeans(posterior)
  # take a sample of 20 posterior draws
  keep <- sample(nrow(posterior), nDraws)
  samp_draws <- posterior[keep, ] 
  
  results <- list(posterios = posterior, posterior_means = post_means, 
                  keep = keep, draws = samp_draws)
  return(results)
  
}