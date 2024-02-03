library(tidyverse)

##### Auxiliary Function to extract fixed effects and credible intervals #####
getFixefQuantREG <- function(fits,
                             probs,
                             robust = TRUE,
                             metric,
                             params = c("Intercept", "Slope"),
                             dimension,
                             Q,
                             QuantLevel, 
                             type) { 
  
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

    fix$Type <- type
    fix$metric <- metric 
    fix$params <- params 
    fix$dimension <- dimension 
    fix$Quantile <- QuantLevel[i] 
    fix$Q <- Q
    fixLst[[i]] <- fix 
    
  } 
  
  fixdt <- do.call(rbind, fixLst)
  rownames(fixdt) <- NULL 
  
  fixdt <- fixdt %>% 
    select(Type, dimension, metric, params, Quantile, Q, everything())
 
  return(fixdt) 
  
}

########### --------- Parameter Estimations  Quantile Regression Distance Based Metrics -------- ###########

##### Phylogenetic dimension ##### 
### Inits
## Probabilities 
Probs <- c(0.025, 0.05, 0.11, 0.25, 0.5, 0.75, 0.89, 0.95, 0.975)
## Quantile level
QuantileLevel <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

### MPD
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_MPD_0595.RData") 

quant_mpd <- getFixefQuantREG(
  fits = quantREG_mpd, 
  probs = Probs, 
  robust = TRUE,
  metric = "MPD",
  dimension = "phylogeny", 
  Q = "Qna",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### MPDz
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_MPDz_0595.RData") 

quant_mpdz <- getFixefQuantREG(
  fits = quantREG_mpdz, 
  probs = Probs, 
  robust = TRUE,
  metric = "MPDz",
  dimension = "phylogeny", 
  Q = "Qna",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### PD
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_PD_0595.RData") 

quant_pd <- getFixefQuantREG(
  fits = quantREG_pd, 
  probs = Probs, 
  robust = TRUE,
  metric = "PD",
  dimension = "phylogeny", 
  Q = "Qna",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### PDz
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_PDz_0595.RData") 

quant_pdz <- getFixefQuantREG(
  fits = quantREG_pdz, 
  probs = Probs, 
  robust = TRUE,
  metric = "PDz",
  dimension = "phylogeny", 
  Q = "Qna",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### qDPM_0
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_qDPM0_0595.RData") 

quant_qdpm0 <- getFixefQuantREG(
  fits = quantREG_qdpm0, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDPM_q0",
  dimension = "phylogeny", 
  Q = "Q0",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### qDPM_1
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_qDPM1_0595.RData") 

quant_qdpm1 <- getFixefQuantREG(
  fits = quantREG_qdpm1, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDPM_q1",
  dimension = "phylogeny", 
  Q = "Q1",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### qDPM_2
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_phylo_qDPM2_0595.RData") 

quant_qdpm2 <- getFixefQuantREG(
  fits = quantREG_qdpm2, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDPM_q2",
  dimension = "phylogeny", 
  Q = "Q2",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

## Combine estimates
quantPHY_estimates_DIS <- bind_rows(quant_mpd, 
                                    quant_mpdz, 
                                    quant_pd, 
                                    quant_pdz, 
                                    quant_qdpm0, 
                                    quant_qdpm1, 
                                    quant_qdpm2)

write_csv(quantPHY_estimates_DIS, 
          file = "Results/Regressions/Reanalyses/FixEffects_table_Phylo_Quantile_results_Distance.csv")

########### Trait dimension ########### 
### Inits
## Probabilities 
Probs <- c(0.025, 0.05, 0.11, 0.25, 0.5, 0.75, 0.89, 0.95, 0.975)
## Quantile level
QuantileLevel <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

### MTD
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_MTD_0595.RData") 

quant_mtd <- getFixefQuantREG(
  fits = quantREG_mtd, 
  probs = Probs, 
  robust = TRUE,
  metric = "MTD",
  dimension = "trait", 
  Q = "Qna",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### MTDz
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_MTDz_0595.RData") 

quant_mtdz <- getFixefQuantREG(
  fits = quantREG_mtdz, 
  probs = Probs, 
  robust = TRUE,
  metric = "MTDz",
  dimension = "trait", 
  Q = "Qna",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### qDTM_0
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_qDTM0_0595.RData") 

quant_qdtm0 <- getFixefQuantREG(
  fits = quantREG_qdtm0, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDTM_q0",
  dimension = "trait", 
  Q = "Q0",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### qDTM_1
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_qDTM1_0595.RData") 

quant_qdtm1 <- getFixefQuantREG(
  fits = quantREG_qdtm1, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDTM_q1",
  dimension = "trait", 
  Q = "Q1",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

### qDTM_2
load("Results/Regressions/Reanalyses/QuantileBLMM/DIS/quant_trait_qDTM2_0595.RData") 

quant_qdtm2 <- getFixefQuantREG(
  fits = quantREG_qdtm2, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDTM_q2",
  dimension = "trait", 
  Q = "Q2",
  QuantLevel = QuantileLevel, 
  type = "Distance"
)

## Combine estimates
quantTRT_estimates_DIS <- bind_rows(quant_mtd, 
                                    quant_mtdz, 
                                    quant_qdtm0, 
                                    quant_qdtm1, 
                                    quant_qdtm2)

write_csv(quantTRT_estimates_DIS, 
          file = "Results/Regressions/Reanalyses/FixEffects_table_Trait_Quantile_results_Distance.csv")

##### Combine all results ###### 

quantPHYLO_estimates_DIS <- read_csv("Results/Regressions/Reanalyses/FixEffects_table_Phylo_Quantile_results_Distance.csv")
quantTRAIT_estimates_DIS <- read_csv("Results/Regressions/Reanalyses/FixEffects_table_Trait_Quantile_results_Distance.csv")

quant_ALL_DIS <- bind_rows(quantPHYLO_estimates_DIS, 
                           quantTRAIT_estimates_DIS)

write_csv(quant_ALL_DIS, 
          file = "Results/Regressions/Reanalyses/FixEffects_table_Quantile_results_DIS.csv")

########### --------- Parameter Estimations  Quantile Regression SS Based Metrics -------- ###########

##### Auxiliary Function to extract fixed effects and credible intervals #####
getFixefQuantREG <- function(fits,
                             probs,
                             robust = TRUE,
                             metric,
                             params = c("Intercept", "Slope"),
                             dimension,
                             Q,
                             QuantLevel, 
                             type) { 
  
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
    
    fix$Type <- type
    fix$metric <- metric 
    fix$params <- params 
    fix$dimension <- dimension 
    fix$Quantile <- QuantLevel[i] 
    fix$Q <- Q
    fixLst[[i]] <- fix 
    
  } 
  
  fixdt <- do.call(rbind, fixLst)
  rownames(fixdt) <- NULL 
  
  fixdt <- fixdt %>% 
    select(Type, dimension, metric, params, Quantile, Q, everything())
  
  return(fixdt) 
  
}

##### Taxonomic dimension #####

### Inits
## Probabilities 
Probs <- c(0.025, 0.05, 0.11, 0.25, 0.5, 0.75, 0.89, 0.95, 0.975)
## Quantile level
QuantileLevel <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

### Richness 
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/quant_taxo_SR_SS_0595.RData")

quant_SR <- getFixefQuantREG( 
  fits = quantREG_sr, 
  probs = Probs, 
  robust = TRUE,
  metric = "S",
  dimension = "taxonomy", 
  Q = "Q0",
  QuantLevel = QuantileLevel, 
  type = "SS"
  )

## Shannon
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/quant_taxo_H_SS_0595.RData")

quant_H <- getFixefQuantREG( 
  fits = quantREG_h, 
  probs = Probs, 
  robust = TRUE,
  metric = "H",
  dimension = "taxonomy", 
  Q = "Q1",
  QuantLevel = QuantileLevel, 
  type = "SS"
)

## Simpson
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/quant_taxo_D_SS_0595.RData")

quant_D <- getFixefQuantREG( 
  fits = quantREG_d, 
  probs = Probs, 
  robust = TRUE,
  metric = "D",
  dimension = "taxonomy", 
  Q = "Q2",
  QuantLevel = QuantileLevel, 
  type = "SS"
)

## Combine estimates
quantTAX_estimates_SS <- bind_rows(quant_SR, 
                                   quant_H, 
                                   quant_D)

write_csv(quantTAX_estimates_SS, 
          file = "Results/Regressions/Reanalyses/FixEffects_table_Taxonomy_Quantile_results_SS.csv")

##### Phylogenetic dimension #####

### Inits
## Probabilities 
Probs <- c(0.025, 0.05, 0.11, 0.25, 0.5, 0.75, 0.89, 0.95, 0.975)
## Quantile level
QuantileLevel <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

### Richness 
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/phylo_q0DPM_SS_quantile.RData")

quant_qdpm0 <- getFixefQuantREG( 
  fits = fit_q0DPM_SS, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDPM_q0",
  dimension = "phylogeny", 
  Q = "Q0",
  QuantLevel = QuantileLevel, 
  type = "SS"
)

## Shannon
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/phylo_q1DPM_SS_quantile.RData")

quant_qdpm1 <- getFixefQuantREG( 
  fits = fit_q1DPM_SS, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDPM_q1",
  dimension = "phylogeny", 
  Q = "Q1",
  QuantLevel = QuantileLevel, 
  type = "SS"
)

## Simpson
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/phylo_q2DPM_SS_quantile.RData")

quant_qdpm2 <- getFixefQuantREG( 
  fits = fit_q2DPM_SS, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDPM_q2",
  dimension = "phylogeny", 
  Q = "Q2",
  QuantLevel = QuantileLevel, 
  type = "SS"
)

## Combine estimates
quantPHYLO_estimates_SS <- bind_rows(quant_qdpm0, 
                                     quant_qdpm1, 
                                     quant_qdpm2)

write_csv(quantPHYLO_estimates_SS, 
          file = "Results/Regressions/Reanalyses/FixEffects_table_Phylogeny_Quantile_results_SS.csv")

##### Trait dimension #####

### Inits
## Probabilities 
Probs <- c(0.025, 0.05, 0.11, 0.25, 0.5, 0.75, 0.89, 0.95, 0.975)
## Quantile level
QuantileLevel <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)

### Richness 
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/trait_q0DTM_SS_quantile.RData")

quant_qdtm0 <- getFixefQuantREG( 
  fits = fit_q0DTM_SS, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDTM_q0",
  dimension = "trait", 
  Q = "Q0",
  QuantLevel = QuantileLevel, 
  type = "SS"
)

## Shannon
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/trait_q1DTM_SS_quantile.RData")

quant_qdtm1 <- getFixefQuantREG( 
  fits = fit_q1DTM_SS, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDTM_q1",
  dimension = "trait", 
  Q = "Q1",
  QuantLevel = QuantileLevel, 
  type = "SS"
)

## Simpson
load("Results/Regressions/Reanalyses/QuantileBLMM/SS/trait_q2DTM_SS_quantile.RData")

quant_qdtm2 <- getFixefQuantREG( 
  fits = fit_q2DTM_SS, 
  probs = Probs, 
  robust = TRUE,
  metric = "qDTM_q2",
  dimension = "trait", 
  Q = "Q2",
  QuantLevel = QuantileLevel, 
  type = "SS"
)

## Combine estimates
quantTRAIT_estimates_SS <- bind_rows(quant_qdtm0, 
                                     quant_qdtm1, 
                                     quant_qdtm2)

write_csv(quantTRAIT_estimates_SS, 
          file = "Results/Regressions/Reanalyses/FixEffects_table_Trait_Quantile_results_SS.csv")

##### Combine all results based on SS #####

quantTAXO_estimates_SS <- read_csv("Results/Regressions/Reanalyses/FixEffects_table_Taxonomy_Quantile_results_SS.csv") 
quantPHYLO_estimates_SS <- read_csv("Results/Regressions/Reanalyses/FixEffects_table_Phylogeny_Quantile_results_SS.csv") 
quantTRAIT_estimates_SS <- read_csv("Results/Regressions/Reanalyses/FixEffects_table_Trait_Quantile_results_SS.csv") 

quant_ALL_SS <- bind_rows(quantTAXO_estimates_SS, 
                          quantPHYLO_estimates_SS, 
                          quantTRAIT_estimates_SS)

write_csv(quant_ALL_SS, 
          file = "Results/Regressions/Reanalyses/FixEffects_table_Quantile_results_SS.csv")
