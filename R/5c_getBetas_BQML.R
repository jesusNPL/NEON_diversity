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

########### --------- Estimations Quantile regression Phylogeny -------- ###########
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_MPDz_site_0595.RData")
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_qDTM0_site_0595.RData")
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_qDTM1_site_0595.RData")
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_qDTM2_site_0595.RData")

quant_mpd <- getFixefQuantREG(
  fits = quantREG_mpd, robust = TRUE,
  metric = "MPDz",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

quant_qdtm0 <- getFixefQuantREG(
  fits = quantREG_qdtm0, robust = TRUE,
  metric = "qDTM_q0",
  dimension = "phylogeny", Q = "Q0",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

quant_qdtm1 <- getFixefQuantREG(
  fits = quantREG_qdtm1, robust = TRUE,
  metric = "qDTM_q1",
  dimension = "phylogeny", Q = "Q1",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

quant_qdtm2 <- getFixefQuantREG(
  fits = quantREG_qdtm2, robust = TRUE,
  metric = "qDTM_q2",
  dimension = "phylogeny", Q = "Q2",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

## Combine estimates
quantPHY_estimates <- rbind(quant_mpd, quant_qdtm0, 
                            quant_qdtm1, quant_qdtm2)

write.csv(quantPHY_estimates, 
          file = "Results/Regressions/phylo-spec/Estimations_FIXEF_Quantphylo_0595.csv")

########### --------- Estimations Quantile regression Trait -------- ###########
load("Results/Regressions/QuantileREG/Trait/quant_trait_MFDz_site_0595.RData")
load("Results/Regressions/QuantileREG/Trait/quant_trait_qDTM0_site_0595.RData")
load("Results/Regressions/QuantileREG/Trait/quant_trait_qDTM1_site_0595.RData")
load("Results/Regressions/QuantileREG/Trait/quant_trait_qDTM2_site_0595.RData")

quant_mfd <- getFixefQuantREG(
  fits = quantREG_mfd, robust = TRUE,
  metric = "MPDz",
  dimension = "trait", Q = "Q13",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

quant_qdtm0 <- getFixefQuantREG(
  fits = quantREG_qdtm0, robust = TRUE,
  metric = "qDTM_q0",
  dimension = "trait", Q = "Q0",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

quant_qdtm1 <- getFixefQuantREG(
  fits = quantREG_qdtm1, robust = TRUE,
  metric = "qDTM_q1",
  dimension = "trait", Q = "Q1",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

quant_qdtm2 <- getFixefQuantREG(
  fits = quantREG_qdtm2, robust = TRUE,
  metric = "qDTM_q2",
  dimension = "trait", Q = "Q2",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

## Combine estimates
quantTRT_estimates <- rbind(quant_mfd, quant_qdtm0, 
                            quant_qdtm1, quant_qdtm2)

write.csv(quantTRT_estimates, 
          file = "Results/Regressions/trait-spec/Estimations_FIXEF_Quanttrait_0595.csv")

########### --------- Estimations Quantile regression Taxonomy -------- ###########
load("Results/Regressions/QuantileREG/Trait/quant_trait_SR_site_0595.RData")
quantREG_srmsd <- quantREG_sr
load("Results/Regressions/QuantileREG/Taxonomy/quant_taxo_SR_site_0595.RData")
load("Results/Regressions/QuantileREG/Taxonomy/quant_taxo_H_site_0595.RData")
load("Results/Regressions/QuantileREG/Taxonomy/quant_taxo_D_site_0595.RData")


quant_sr1 <- getFixefQuantREG(
  fits = quantREG_srmsd, robust = TRUE,
  metric = "S1",
  dimension = "taxonomy", Q = "Q13",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

quant_sr2 <- getFixefQuantREG(
  fits = quantREG_sr, robust = TRUE,
  metric = "S2",
  dimension = "taxonomy", Q = "Q0",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

quant_h <- getFixefQuantREG(
  fits = quantREG_h, robust = TRUE,
  metric = "Shannon",
  dimension = "taxonomy", Q = "Q1",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

quant_d <- getFixefQuantREG(
  fits = quantREG_d, robust = TRUE,
  metric = "Simpson",
  dimension = "taxonomy", Q = "Q2",
  QuantLevel = c(0.05, 0.1, 0.90, 0.95)
)

## Combine estimates
quantTRT_estimates <- rbind(quant_sr1, quant_sr2, 
                            quant_h, quant_d)

write.csv(quantTRT_estimates, 
          file = "Results/Regressions/taxo-spec/Estimations_FIXEF_Quanttaxonomy_0595.csv")

