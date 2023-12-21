library(tidyverse)
library(spectrolab)
library(raster)

########## Step 1 - mask spectra data ##########

##### Auxiliary functions #####

### calculate NDVI
NDVI <- function(NIR, RED) {
  (NIR - RED)/(NIR + RED)
}

### NIR mask
NIRmask <- function(NIR_bands) { 
  
  nir <- calc(NIR_bands, mean)
  
  scale_factor <- 10000 
  
  Mask <- nir / scale_factor 
  
  ## Code extracted from https://github.com/akamoske/hypRspec/blob/master/R/brightnessmask.R
  # lets calculate the mean reflectance value for this matrix
  r.quantiles <- quantile(Mask, probs = seq(0, 1, 0.25), na.rm = TRUE)
  
  r.25 <- r.quantiles[[2]] 
  
  r.75 <- r.quantiles[[4]]
  
  r.iqr <- r.75 - r.25
  
  r.bottom.cut <- r.25 - (1.5 * r.iqr)
  
  Mask[Mask <= r.bottom.cut] <- NA 
  
  return(Mask)
  
}

##### Mask HSI using NDVI and NIR #####

### Inits 

## Sites file names
fileNames <- list.files("DATA/Spectra/", pattern = ".csv") 

## Sites names
siteNames <- fileNames %>% 
  tibble(sites = .) %>% 
  separate(sites, c("siteID")) %>% 
  as.data.frame()
  
## Number of sites
nSites <- length(fileNames)

### Start process
for(i in 1:nSites) { 
  
  print(paste0("Loading ", siteNames$siteID[i])) 
  
  ## Load data
  site <- read_csv(paste0("DATA/Spectra/", fileNames[i])) %>% 
    drop_na(plotID)
  
  ## Plots names
  plotNames <- unique(site$plotID)
  ## number of plots
  nPlots <- length(plotNames)
  
  ## Store masked plots
  plot_masked <- list()
  
  for(j in 1:nPlots) { 
    
    print(paste0("Cleanning and/or masking HSI data for plot ", plotNames[j]))  
    
    ### Filter plots 
    pt <- site %>% 
      dplyr::filter(plotID == plotNames[j]) 
    
    ### Create raster from plot 
    pt_ras <- raster::brick(sp::SpatialPixelsDataFrame(points = pt[, 1:2], pt)) 
    
    ##### Calculate NDVI #####
    # select bands to use in calculation (red, NIR) - bands c(58, 90) in full NEON hyperspectral dataset 
    # https://www.neonscience.org/resources/learning-hub/tutorials/create-raster-stack-hsi-hdf5-r
    red <- pt_ras[[62]] # red band
    
    nir <- pt_ras[[94]] # nir band
    
    ndvi <- NDVI(NIR = nir, RED = red)
    
    ## global NVDI threshold
    ndvi_threshold <- ndvi 
    
    ndvi_threshold[ndvi_threshold <= 0.5] <- NA
    
    ## threshold for plots with low vegetation coverage
    ndvi_threshold_2 <- ndvi 
    
    ndvi_threshold_2[ndvi_threshold_2 <= 0.2] <- NA
    
    ##### Mask plot using NDVI ##### 
    pt_test_ndvi <- mask(subset(pt_ras, 5:nlayers(pt_ras)), ndvi_threshold) # select only rasters with spectral information
    
    pt_test_ndvi_df <- as.data.frame(pt_test_ndvi) %>% 
      drop_na()
    
    if (nrow(pt_test_ndvi_df) > 200) { # if the plot has more than 200 pixels continue
      
      pt_ras_ndvi <- mask(subset(pt_ras, 5:nlayers(pt_ras)), ndvi_threshold) # select only rasters with spectral information
      
    } else { # else, use threshold for plots with low vegetation coverage
      
      pt_ras_ndvi <- mask(subset(pt_ras, 5:nlayers(pt_ras)), ndvi_threshold_2) # select only rasters with spectral information
      
    }
    
    ### Obtain threshold for NIR
    nir_threshold <- NIRmask(NIR_bands = subset(pt_ras_ndvi, 75:134)) # only NIR bands
    
    ##### Mask plot using NIR #####
    pt_ras_ndvi_nir <- mask(pt_ras_ndvi, nir_threshold) 
    
    pt_df <- as.data.frame(pt_ras_ndvi_nir)
    
    pt_df <- cbind(pt[, c(1, 2, 4)], pt_df)
    
    plot_masked[[j]] <- pt_df
    
  } # end for j
  
  plot_masked <- do.call(rbind, plot_masked)
  
  write_csv(plot_masked, 
            file = paste0("DATA/Spectra/masked/mask_", fileNames[i]))
  
  print(paste0("Masking for site ", siteNames$siteID[i], ", done!"))
  
} # end for i

### End process

########## Step 2 - spectra normalization ##########
library(tidyverse)
library(spectrolab)

### Inits 
## Sites files names
fileNames <- list.files("DATA/Spectra/masked/", pattern = ".csv") 

## Sites names
siteNames <- fileNames %>% 
  tibble(sites = .) %>% 
  separate(sites, c("type", "siteID", "dato")) %>% 
  dplyr::select(siteID) %>% 
  as.data.frame()

## Number of sites
nSites <- length(fileNames)

## Start process
for(i in 1:nSites) { 
  
  print(paste0("Loading site ", siteNames$siteID[i])) 
  
  ## Load data
  site <- read_csv(paste0("DATA/Spectra/masked/", fileNames[i])) %>% 
    drop_na(plotID)
  
  ## Plots
  plotNames <- unique(site$plotID)
  nPlots <- length(plotNames)
  
  ## Store masked plots
  plot_normalization <- list() 
  
  for(j in 1:nPlots) { 
    
    print(paste0("spectra normalization for ", plotNames[j]))
    
    ## Filter plots 
    pt <- site %>% 
      dplyr::filter(plotID == plotNames[j]) 
    
    ### Remove water bands 
    ## long format
    pt_long <- pt %>% 
      pivot_longer(!c(X, Y, plotID), names_to = "band", values_to = "reflectance") %>% 
      separate(band, c("index", "wavelength")) %>%
      mutate(wavelength = parse_number(wavelength)) %>%
      filter(
        !between(wavelength, 1340, 1445), # filter water bands
        !between(wavelength, 1790, 1955)
      ) %>% 
      filter(
        wavelength >= 400, # remove bands with noise
        wavelength <= 2450) 
    
    ### Normalization 
    ## wide format
    pt_wide <- pt_long %>% 
      pivot_wider(id_cols = c(X, Y, plotID), 
                  names_from = wavelength, 
                  values_from = reflectance) 
    
    pt_spec <- spectrolab::as_spectra(pt_wide[, 3:ncol(pt_wide)], name_idx = 1)
    
    ### Normalize spectra
    pt_norm <- normalize(pt_spec)
    
    pt_normalized <- pt_norm %>% 
      as.data.frame() %>% 
      dplyr::select(-c(normalization_magnitude, sample_name)) 
      
    pt_normalized <- cbind(pt_wide[, 1:3], pt_normalized)
    
    plot_normalization[[j]] <- pt_normalized
    
    } # end for j
  
  plot_normalization <- do.call(rbind, plot_normalization)
  
  write_csv(plot_normalization, 
            file = paste0("DATA/Spectra/normalized/normalized_", fileNames[i]))
  
  print(paste0("Masking for site ", siteNames$siteID[i], ", done!"))
  
} # end for i

### End process
