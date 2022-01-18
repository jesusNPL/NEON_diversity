library(tidyverse)
library(rgdal)
library(sp)
library(raster)
library(spatialEco)

### Merge domains and biomes
NEON_domains <- readOGR(dsn = "DATA/NEONDomains_0", 
                        layer = "NEON_Domains")

NEON_domains@data <- NEON_domains@data %>% 
  mutate(dID = c("D01", "D01", "D02", "D03", "D04", "D04", "D05", "D06", "D07", 
                 "D08", "D09", "D10", "D11", "D12", "D13", "D14", "D15", "D16", 
                 "D17", "D18", "D19", "D20"))

eco2017 <- readOGR("../../Downloads/Ecoregions2017/Ecoregions2017.shp")
view(eco2017@data)

biomesNA <- crop(eco2017, NEON_domains)
view(biomesNA)

NEON_biomes <- raster::intersect(NEON_domains, biomesNA)
view(NEON_biomes)

NEON_domains_biomes <- NEON_biomes@data %>% 
  dplyr::select(DomainID, dID, DomainName, ECO_NAME, BIOME_NUM, 
         BIOME_NAME, ECO_BIOME_)

save(NEON_biomes, biomesNA, file = "DATA/Spatial/NEON_TEW_merge.RData")

load("DATA/Spatial/NEON_TEW_merge.RData")

write.csv(NEON_domains_biomes, "DATA/NEON_TEW_combined.csv")

##### Load and combine data #####

NEON_TEW <- read.csv("DATA/NEON_TEW_combined.csv")

NEON_TEW <- NEON_TEW %>% 
  mutate(codeRep = paste0(dID, "_", BIOME_NUM))

NEON_TEW <- NEON_TEW[!duplicated(NEON_TEW$codeRep), ]

NEON_TEW[duplicated(NEON_TEW$codeRep), ]

head(NEON_TEW)

## Load neon plots
NEON_plots <- readOGR(dsn = "DATA/spatial/All_NEON_TOS_Plots_v8", 
                      layer = "All_NEON_TOS_Plot_Polygons_V8")

neon_mt <- NEON_plots@data %>% 
  dplyr::select(plotIDNEON = plotID, domain, domainID, siteID, 
                latitude, longitude, elevation, nlcdClass, soilOrder)

neon_mt <- neon_mt[!duplicated(neon_mt$plotIDNEON), ]

neon_mt[duplicated(neon_mt$plotIDNEON), ]

head(neon_mt)

## combine neon plots with neon maps and biome info
neon_mt <- left_join(neon_mt, NEON_TEW, 
                    by = c("domainID" = "dID")
          )

neon_mt <- neon_mt[!duplicated(neon_mt$plotIDNEON), ]

neon_mt[duplicated(neon_mt$plotIDNEON), ]

head(neon_mt)

## filter neon sites available in our dataset
load("Results/RegDATA/phylo_spec_DIS_NEON.RData")

NEON_sites <- unique(phylo_spec_NEON_table_q0$Site)

neon_mt <- neon_mt[neon_mt$siteID %in% NEON_sites, ]

neon_mt[duplicated(neon_mt$plotIDNEON), ]

##### Prepare data to match with metrics #####
neon_mt <- neon_mt %>% 
  mutate(plots = plotIDNEON) %>% 
  separate(plots, into = c("site", "plot"), sep = "_")

neon_mt <- neon_mt %>% 
  mutate(plotID = str_remove(neon_mt$plot, "[00]")) %>% 
  mutate(plotID2 = paste0("P_", plotID)) %>% 
  mutate(plotID = str_replace(plotID2, "_0", "_")) 

neon_mt <- neon_mt %>% 
  dplyr::select(domainID, domain, siteID, site, 
                plot, plotIDNEON, plotID, plotID2, 
                everything() 
                )
head(neon_mt)

neon_mt[duplicated(neon_mt$plotIDNEON), ]

write.csv(neon_mt, "Results/RegDATA/NEON_metadata.csv")

neon_mt2 <- neon_mt[neon_mt$plotID %in% unique(phylo_spec_NEON_table_q0$plotID), ]

write.csv(neon_mt2, "Results/RegDATA/NEON_metadata_MATCH_plotID.csv")


xx <- left_join(phylo_spec_NEON_table_q0, neon_mt2, 
                by = c("Site" = "site", "plotID"))

yy <- xx %>% 
  filter(nlcdClass == "evergreenForest")

zz <- xx %>% 
  filter(nlcdClass == "grasslandHerbaceous")

ww <- xx %>% 
  filter(nlcdClass == "shrubScrub")

uu <- xx %>% 
  filter(nlcdClass == "mixedForest") 

vv <- xx %>% 
  filter(nlcdClass == "deciduousForest")

yy[duplicated(yy$plotIDNEON), ]

library(brms)

my <- brm(
  bf(qDTM_p ~ qDTM_s + (1 | Site)),
  data = yy,
  control = list(adapt_delta = 0.9),
  chains = 4, iter = 2000, warmup = 1000, cores = 4, 
  backend = "cmdstanr"
)

mz <- brm(
  bf(qDTM_p ~ qDTM_s), # + (1 | Site)),
  data = zz,
  control = list(adapt_delta = 0.9),
  chains = 4, iter = 2000, warmup = 1000, cores = 4,
  backend = "cmdstanr"
)

mw <- brm(
  bf(qDTM_p ~ qDTM_s + (1 | Site)),
  data = ww,
  control = list(adapt_delta = 0.9),
  chains = 4, iter = 2000, warmup = 1000, cores = 4, 
  backend = "cmdstanr"
)

mu <- brm(
  bf(qDTM_p ~ qDTM_s + (1 | Site)),
  data = uu,
  control = list(adapt_delta = 0.9),
  chains = 4, iter = 2000, warmup = 1000, cores = 4, 
  backend = "cmdstanr"
)

mv <- brm(
  bf(qDTM_p ~ qDTM_s), # + (1 | Site)),
  data = vv,
  control = list(adapt_delta = 0.9),
  chains = 4, iter = 2000, warmup = 1000, cores = 4, 
  backend = "cmdstanr"
)

library("bayesplot")
theme_set(bayesplot::theme_default(base_family = "sans"))

posterior <- posterior_samples(my, "^b")
post_means <- colMeans(posterior)
# take a sample of 20 posterior draws
keep <- sample(nrow(posterior), 50)
samp_draws <- posterior[keep, ]

yy %>%
  ggplot(aes(x = qDTM_s, y = qDTM_p)) +
  geom_abline( # 20 posterior draws of the regression line
    intercept = samp_draws[, 1],
    slope = samp_draws[, 2],
    color = "#9497eb",
    size = 1
  ) +
  geom_abline( # posterior mean regression line
    intercept = post_means[1],
    slope = post_means[2],
    color = "#1c35c4",
    size = 1
  ) +
  geom_point(alpha = 0 / 10
  ) +
  # ylim(range(d_pd$PD)) +
  ylim(0, 20) +
  # xlim(1, 20) +
  labs(x = "qDTM - Spectra", y = "qDTM - Phylogeny") +
  labs(title = bquote(R^2 ~ "=" ~ "0.45[0.41:0.48]"))


neon_metadata <- neon_mt %>% 
  dplyr::select(plotIDNEON, domain, domainID, siteID, elevation) %>% 
  group_by(domain, domainID, siteID) %>% 
  summarize(elev = mean(elevation))

NEON_TEW_sel <- NEON_TEW[NEON_TEW$dID %in% neon_metadata$domainID, ]

NEON_TEW_sel <- NEON_TEW_sel %>% 
  dplyr::select(DomainID, dID, DomainName, BIOME_NUM, BIOME_NAME) %>% 
  mutate(BIOME_NUM = paste0("B_", BIOME_NUM), domainID = dID) %>% 
  group_by(BIOME_NUM, BIOME_NAME, DomainName, domainID) %>% 
  summarize(val = mean(DomainID)) 

NEON_META <- full_join(neon_metadata, NEON_TEW_sel, by = "domainID")

NEON_META$habitat <- case_when(NEON_META$BIOME_NAME == "Temperate Broadleaf & Mixed Forests" ~ "forest", 
          NEON_META$BIOME_NAME == "Temperate Grasslands, Savannas & Shrublands" ~ "open", 
          NEON_META$BIOME_NAME == "Tropical & Subtropical Moist Broadleaf Forests" ~ "forest", 
          NEON_META$BIOME_NAME == "Mangroves" ~ "other", 
          NEON_META$BIOME_NAME == "Tropical & Subtropical Dry Broadleaf Forests" ~ "forest", 
          NEON_META$BIOME_NAME == "Flooded Grasslands & Savannas" ~ "open", 
          NEON_META$BIOME_NAME == "Deserts & Xeric Shrublands" ~ "open", 
          NEON_META$BIOME_NAME == "Temperate Conifer Forests" ~ "forest", 
          NEON_META$BIOME_NAME == "Mediterranean Forests, Woodlands & Scrub" ~ "forest", 
          NEON_META$BIOME_NAME == "Tropical & Subtropical Coniferous Forests" ~ "forest", 
          NEON_META$BIOME_NAME == "N/A" ~ "tundra", 
          NEON_META$BIOME_NAME == "Tundra" ~ "tundra", 
          NEON_META$BIOME_NAME == "Boreal Forests/Taiga" ~ "forest", 
          NEON_META$BIOME_NAME == "Tropical & Subtropical Grasslands, Savannas & Shrublands" ~ "open"
          )

view(NEON_META)

write.csv(NEON_META, "Results/RegDATA/NEON_DOMAINS_BIOMES_HABITAT.csv")

NEON_META <- read.csv("Results/RegDATA/NEON_DOMAINS_BIOMES_HABITAT.csv")

NEON_META %>% 
  group_by(siteID, habitat) %>% 
  summarize(ele = mean(elev)) %>% 
  view()

?conditional_effects
