makeTable_taxonomy <- function(div1, div2) {
  # Vector used to rename the taxonomic diversities
  name1 <- c("SrObs", "SrObs_SErr", "ShannonObs", "ShannonObs_SErr", 
             "SimpsonObs", "SimpsonObs_SErr", 'H_Obs', "SimpObs", "invSimpObs", 
             "unbiasSimpObs", "FisherObs", "S_Obs", "J_Obs", "plotID") 
  
  name2 <- c("SrSpec", "SrSpec_SErr", "ShannonSpec", "ShannonSpec_SErr", 
             "SimpsonSpec", "SimpsonSpec_SErr", 'H_Spec', "SimpSpec", "invSimpSpec", 
             "unbiasSimpSpec", "FisherSpec", "S_Spec", "J_Spec", "plotID") 
  # Rename columns
  names(div1) <- name1
  names(div2) <- name2
  
  tab <- full_join(div1, div2, by = "plotID") 
  tab <- tab[, c(14, 1:13, 15:27)]
  
  obs <- tab[, c(1, 2, 4, 6, 8:14, 15, 17, 19, 21:27)]
  SErr <- tab[, c(1, 3, 5, 7, 16, 18, 20)]
  
  res <- list(divOBS = obs, divSErr = SErr)
  return(res)
}
