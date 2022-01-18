

harv <- read.csv("Dropbox/Macrosystems_NEON/DATA/Harvard_NEON/Data/HARV_spectra_grid.csv")

harvPlots <- unique(harv$plotID)[1:33]
nPlots <- length(harvPlots)

## make rasters
harv_ras <- demon_SPEC_to_RASTER(spectra = harv, 
                                 plotNames = harvPlots, 
                                 nPlots = nPlots)

## calculate SAM
harv_sam <- distSAM(spectraRAS = harv_ras, 
                    plotNames = harvPlots, 
                    nPlots = nPlots)

### make metrics
harv_metrics <- demon_DivDistance_SAM(samList = harv_sam, 
                                 Q = 0, nPlots = nPlots, 
                                 plotNames = harvPlots
)

