# Make community data

make_SSppComm <- function(spectra) {

  samp <- spectra %>% 
    count(plotID, specSpp)
  samp <- samp[, c(1, 3, 2)]
  comm <- picante::sample2matrix(samp)
  
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
  
  comm <- round(comm + 0.5)
  
  alphas$H <- diversity(comm)
  alphas$Simp <- diversity(comm, "simpson")
  alphas$invSimp <- diversity(comm, "inv")
  
  if (!identical(all.equal(comm, round(comm)), TRUE)) {
    comm2 <- decostand(comm, "max")
    comm2 <- wisconsin(comm2)
  
    ## Fisher alpha
    alphas$Fisher_alpha <- fisher.alpha(comm2)
    ## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
    alphas$unbiasSimp <- rarefy(comm2, 2) - 1
    
  } else { 
    ## Fisher alpha 
    comm2 <- decostand(comm, "total")
    comm2 <- wisconsin(comm2)
    ## Fisher alpha
    alphas$Fisher_alpha <- fisher.alpha(comm2)
    ## Unbiased Simpson (Hurlbert 1971, eq. 5) with rarefy:
    alphas$unbiasSimp <- rarefy(comm2, 2) - 1
    } 
  ## Species richness (S) and Pielou's evenness (J):
  alphas$S <- specnumber(comm) ## rowSums(BCI > 0) does the same...
  alphas$J_Pielou <- alphas$H/log(alphas$S)
  ## Plot all
  alphas$plotID <- plotNames
  alphas$Site <- site
  
  print("Alpha diversity computed at plot level... ")
  
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

