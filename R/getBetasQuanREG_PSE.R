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
### PSR
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_PSR_MSD_0595.RData")
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_PSR_SD_0595.RData")
### SR
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_SR_MSD_0595.RData")
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_SR_SD_0595.RData")
### PSE
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_PSE_MSD_0595.RData")
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_PSE_SD_0595.RData")
### PSV
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_PSV_MSD_0595.RData")
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_PSV_SD_0595.RData")

### PD 
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_PD_0595.RData")
load("Results/Regressions/QuantileREG/Phylo/quant_phylo_PDz_0595.RData")

##### MSD
quant_psr_msd <- getFixefQuantREG(
  fits = quantREG_msd_psr, robust = TRUE,
  metric = "PSR",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.1, 0.75, 0.90, 0.95)
)

quant_sr_msd <- getFixefQuantREG(
  fits = quantREG_msd_sr, robust = TRUE,
  metric = "SR",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.1, 0.75, 0.90, 0.95)
)

quant_pse_msd <- getFixefQuantREG(
  fits = quantREG_msd_pse, robust = TRUE,
  metric = "PSE",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.1, 0.75, 0.90, 0.95)
)

quant_psv_msd <- getFixefQuantREG(
  fits = quantREG_msd_psv, robust = TRUE,
  metric = "PSV",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.1, 0.75, 0.90, 0.95)
)

quant_msd_pse <- rbind(quant_pse_msd, quant_psv_msd, 
                       quant_psr_msd, quant_sr_msd) 

write_csv(quant_msd_pse, 
          file = "Results/Regressions/phylo-spec/Estimations_FIXEF_Quant_PSE_MSD.csv")

##### PD

quant_psr_sd <- getFixefQuantREG(
  fits = quantREG_sd_psr, robust = TRUE,
  metric = "PSR_SD",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.1, 0.75, 0.90, 0.95)
)

quant_sr_sd <- getFixefQuantREG(
  fits = quantREG_sd_sr, robust = TRUE,
  metric = "SR_SD",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.1, 0.75, 0.90, 0.95)
)

quant_pse_sd <- getFixefQuantREG(
  fits = quantREG_sd_pse, robust = TRUE,
  metric = "PSE_SD",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.1, 0.75, 0.90, 0.95)
)

quant_psv_sd <- getFixefQuantREG(
  fits = quantREG_sd_psv, robust = TRUE,
  metric = "PSV_SD",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.1, 0.75, 0.90, 0.95)
)

quant_pd_sd <- getFixefQuantREG(
  fits = quantREG_pd, robust = TRUE,
  metric = "PD_SD",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
)

quant_pdz_sd <- getFixefQuantREG(
  fits = quantREG_pdz, robust = TRUE,
  metric = "PDz_SD",
  dimension = "phylogeny", Q = "Q13",
  QuantLevel = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
)

quant_sd_pse <- rbind(quant_pse_sd, quant_psv_sd, 
                       quant_psr_sd, quant_sr_sd, 
                      quant_pd_sd, quant_pdz_sd) 

write_csv(quant_sd_pse, 
          file = "Results/Regressions/phylo-spec/Estimations_FIXEF_Quant_PSE_SD.csv")

