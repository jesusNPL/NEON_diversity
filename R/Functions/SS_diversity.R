##### Function to make community data based on spectral species #####

make_SSppComm <- function(SS_spectra) {

  ### Select and count SS per plot
  samp <- SS_spectra %>% 
    count(plotID, SS) %>% 
    select(plotID, n, SS)

  ### create community composition
  # Using pivot wider
  comm <- samp %>% 
    pivot_wider(id_cols = plotID, 
                names_from = SS, 
                values_from = n, 
                values_fill = 0) 
  
  # Using picante package
  #comm <-  samp %>% 
    #picante::sample2matrix(.)
  
  return(comm) 
  
}


##### Calculate SS diversity #####
# comm = community data matrix
# q1 to q3 = different levels of q

demon_SSppDIV <- function(comm, q1, q2, q3, plotNames, nPlots, site) {
  
  alphas <- data.frame(matrix(ncol = 6, nrow = nPlots))
  
  for(i in 1:nPlots) { 
    
    # Q = 0
    alfa1 <- d(comm[i, ], q = q1, boot = TRUE, boot.arg = list(num.iter = 1000)) 
    
    alphas[i, 1] <- alfa1[[1]]
    alphas[i, 2] <- alfa1[[2]]
    
    # Q = 1
    alfa2 <- d(comm[i, ], q = q2, boot = TRUE, boot.arg = list(num.iter = 1000)) 
    
    alphas[i, 3] <- alfa2[[1]]
    alphas[i, 4] <- alfa2[[2]]
    
    # Q = 2
    alfa3 <- d(comm[i, ], q = q3, boot = TRUE, boot.arg = list(num.iter = 1000)) 
    
    alphas[i, 5] <- alfa3[[1]]
    alphas[i, 6] <- alfa3[[2]]
    
  } 
  
  names(alphas) <- c("SRich", "SRichSErr", "Shannon", "ShannonSErr", "Simpson", "SimpsonSErr")
  
  ## Add site and plot names
  alphas$Site <- site 
  alphas$plotID <- plotNames 
  
  print("Alpha diversity computed at plot level... ")
  
  ##### Site scale calculations #####
  ### Inits 
  alphaSite <- data.frame(matrix(ncol = 2, nrow = 3))
  betaSite <- data.frame(matrix(ncol = 2, nrow = 3))
  gammaSite <- data.frame(matrix(ncol = 2, nrow = 3))
  
  diver <- c("alpha", "beta", "gamma") 
  
  qs <- c(0, 1, 2)
  
  for(j in 1:length(qs)) {
    
    # Alpha
    alpha <- d(comm, lev = diver[1], q = qs[j], boot = TRUE, 
               boot.arg = list(num.iter = 1000)) 
    
    alphaSite[j, 1] <- alpha[[1]] 
    alphaSite[j, 2] <- alpha[[2]]
    
    # Beta
    beta <- d(comm, lev = diver[2], q = qs[j], boot = TRUE, 
               boot.arg = list(num.iter = 1000)) 
    
    betaSite[j, 1] <- beta[[1]] 
    betaSite[j, 2] <- beta[[2]]
    
    # Gamma
    gamma <- d(comm, lev = diver[3], q = qs[j], boot = TRUE, 
               boot.arg = list(num.iter = 1000)) 
    
    gammaSite[j, 1] <- gamma[[1]] 
    gammaSite[j, 2] <- gamma[[2]]
    
  } 
  
  print("Alpha-Beta-Gamma diversities computed at site level")
  
  # Alpha
  alphaSite$q <- c("q0_Richness", "q1_Shannon", "q2_Simpson")
  alphaSite$Div <- "Alpha"
  names(alphaSite) <- c("D.Value", "StdErr", "q", "Diversity")
  # Beta
  betaSite$q <- c("q0_Richness", "q1_Shannon", "q2_Simpson")
  betaSite$Div <- "Beta"
  names(betaSite) <- c("D.Value", "StdErr", "q", "Diversity")
  # Gamma
  gammaSite$q <- c("q0_Richness", "q1_Shannon", "q2_Simpson")
  gammaSite$Div <- "Gamma"
  names(gammaSite) <- c("D.Value", "StdErr", "q", "Diversity")
  
  sitediv <- rbind(alphaSite, betaSite, gammaSite)
  sitediv$Site <- site
  diversities <- list("commDiv" = alphas, "siteDiv" = sitediv)
  
  return(diversities)
  
}

