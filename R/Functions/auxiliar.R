##### Combine tables taxonomy #####

makeTable_taxonomy <- function(div1, div2) {
  # Vectors used to rename the taxonomic diversities
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


##### Combine table phylogeny #####

makeTable_phylogeny <- function(div1, div2) {
  # Vectors used to rename the taxonomic diversities
  name1 <- c("plotID", "PD", "PDz", "MPD", "MPDz", "Rao_Phylo", "SR", 
             "qHill", "M_Phylo", "mPrime_Phylo", "qHt_Phylo", 
             "qEt_Phylo" , "qDT_Phylo", "qDTM_Phylo" ) 
  
  name2 <- c("plotID", "SD", "MSD", "MSm", "MNSD", "MxNSD", "SNTD", "Range", "nPixels", 
             "qHill", "M_Spec", "mPrime_Spec" , "qHt_Spec", 
             "qEt_Spec", "qDT_Spec", "qDTM_Spec", 
             "SDivergence", "SEvenness", "Rao_Spec") 
  # Rename columns
  names(div1) <- name1
  names(div2) <- name2
  
  # Select columns
  div1 <- div1 %>% 
    select(plotID, PD, PDz, MPD, MPDz, Rao_Phylo, qHill, 
           M_Phylo, mPrime_Phylo, qHt_Phylo, qEt_Phylo, qDT_Phylo, qDTM_Phylo)
  div2 <- div2 %>% 
    select(plotID, SD, MSD, MSm, Rao_Spec, Range, qHill, 
           M_Spec, mPrime_Spec, qHt_Spec, qEt_Spec, qDT_Spec, qDTM_Spec)
  
  tab <- full_join(div1, div2, by = c("plotID", "qHill")) 
  return(tab)
}
