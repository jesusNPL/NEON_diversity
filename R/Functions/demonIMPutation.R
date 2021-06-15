##### Functions to perform phylogenetic imputation #####


makePVR <- function(phylo, multiPhylo) {
  print("Starting phylogenetic decomposition for best tree ...")
  pvr <- PVR::PVRdecomp(phylo)
  
  print("Starting phylogenetic decomposition for sample of trees ...")
  
  pvrx <- pbapply::pblapply(multiPhylo, PVR::PVRdecomp)
  
  PVRs <- list(pvrBest = pvr, pvrSample = pvrx)
}


demon_ImpTrait <- function(traits, traitNames, phylo, multiPhylo, 
                           iters = 100, ntrees = 250, pvr, pvrSamp) { 
  # Select number of trees
  #multiPhylo <- sample(multiPhylo, nPhy)
  
  # Sort traits according tip labels
  trait <- traits[match(phylo$tip.label, rownames(traits)), ]
  
  ## Imputation
  # Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
  pvr <- pvr
  # Extract the PVRs
  pvrs <- pvr@Eigen$vectors
  # Combine traits and PVRs
  traits_pvrs <- cbind(trait, pvrs)
  
  ###### Perform imputation
  # Imputation using missForest (note that this function have other arguments, see details in Stekhoven and Buehlmann 2012)
  imp <- missForest::missForest(traits_pvrs, maxiter = iters, ntree = ntrees, 
                                   variablewise = FALSE)
  
  impTraits <- imp$ximp[, traitNames]
  # OOBerror
  err <- imp$OOBerror
  
  print("Imputation for best tree, done...")
  
  #bestImp <- list(impTraits, err)
  
  ##### Perform imputation using a multiphylo
  
  #pvrx <- list()
  # Extract PVRs
  #for(i in 1:length(multiPhylo)) {
   # print(paste0("Starting phylogenetic decomposition for tree ", i, " ..."))
    #phy <- multiPhylo[[i]] 
    # Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
    #pvrx[[i]] <- PVR::PVRdecomp(phy)
  #}
  
  # Imputed traits
  impTraits_mtphylo <- list()
  # Error associated
  errors <- numeric(length = length(multiPhylo))
  
  # Sort traits according tip labels
  for(j in 1:length(multiPhylo)){
    print(paste0("Initializing imputation using tree ", j, "..."))
    
    phy <- multiPhylo[[j]]
    trait <- traits[match(phy$tip.label, rownames(traits)), ]
    ## Imputation
   
    # Extract the PVRs
    pvrs <- pvrSamp[[j]]@Eigen$vectors
    # Combine traits and PVRs
    traits_pvrs <- cbind(trait, pvrs)
    
    ###### Perform imputation
    # Imputation using missForest (note that this function have other arguments, see details in Stekhoven and Buehlmann 2012)
    imp <- missForest::missForest(traits_pvrs, maxiter = iters, ntree = ntrees, 
                                  variablewise = FALSE)
    
    impTraits_mtphylo[[j]] <- imp$ximp[, traitNames]
    # OOBerror
    errors[j] <- imp$OOBerror
    
    print(paste0("Imputation for tree ", j, " ...done!"))
    
  }
  
  print("Imputation for sample of trees, done!!!")
  
  results <- list(impTraits = impTraits, NRMSE = err, 
                  impMulti = impTraits_mtphylo, NRMSEs = errors)
  return(results)
  
}
