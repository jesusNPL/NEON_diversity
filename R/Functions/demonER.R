# Function to extract hypothesis 
demonER <- function(fit, metric, alpha, quantile) {
  
  dtmat <- data.frame(brms::fixef(fit, robust = TRUE))
  
  slope <- dtmat["Estimate"][2, ]
  param <- rownames(dtmat)[2]
  
  if (alpha == 0.1) { 
    
    if (slope > 0) {
      hyp <- paste0(param, " > 0")
      h <- brms::hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    } else {
      hyp <- paste0(param, " < 0")
      h <- brms::hypothesis(fit, hyp, class = "b", alpha = alpha)$hypothesis
    }
  } else {
    
    if (slope > 0) {
      hyp <- paste0(param, " > 0")
      h <- brms::hypothesis(fit, hyp, class = "b", alpha = 0.05)$hypothesis
    } else {
      hyp <- paste0(param, " < 0")
      h <- brms::hypothesis(fit, hyp, class = "b", alpha = 0.05)$hypothesis
    } 
  }
  
  dtmat$Parameters <- c("Intercept", "Slope")
  dtmat$ER <- h$Evid.Ratio
  dtmat$PostProb <- h$Post.Prob
  dtmat$Star <- h$Star 
  dtmat$Alpha <- alpha
  dtmat$Metric <- metric 
  dtmat$Quantile <- quantile
  
  rownames(dtmat) <- NULL
  
  return(dtmat)
}
