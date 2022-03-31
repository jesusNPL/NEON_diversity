library(tidyverse)

source("R/NEON_diversity/R/Functions/demonRarefaction.R")

##### Load data #####
load("DATA/matchTraitComm/matchTraitComm.RData") 
sites <- names(matchedTrait$commNEON)

##### Run rarefaction for q0, q1, q2 #####
sampleCompletenessNEON <- demoniNEXT(samples = matchedTrait$commNEON, 
                                     Q = c(0, 1, 2), 
                                     nBoots = 1000)

save(sampleCompletenessNEON, 
     file = "output/NEON_sample_completeness.RData")

##### Run rarefaction for q0 #####
sampleCompletenessNEON_q0 <- demoniNEXT(samples = matchedTrait$commNEON, 
                                        Q = 0,  
                                        nBoots = 1000)
save(sampleCompletenessNEON_q0, 
     file = "output/NEON_sample_completeness_q0.RData")

##### extract Estimators ##### 
load("output/NEON_sample_completeness.RData")

## SR
esti_SR_NEON <- sampleCompletenessNEON$Estimators %>% 
  filter(Diversity == "Species richness") 

view(esti_SR_NEON)

cor.test(esti_SR_NEON$Observed, esti_SR_NEON$Estimator)

esti_SR_NEON %>% 
  filter(Estimator < 80) %>% 
  ggplot(aes(x = Observed, y = Estimator)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0)

## Shannon
esti_H_NEON <- sampleCompletenessNEON$Estimators %>% 
  filter(Diversity == "Shannon diversity")

cor.test(esti_H_NEON$Observed, esti_H_NEON$Estimator)
plot(esti_H_NEON$Observed, esti_H_NEON$Estimator)

## Simpson
esti_D_NEON <- sampleCompletenessNEON$Estimators %>% 
  filter(Diversity == "Simpson diversity")

cor.test(esti_D_NEON$Observed, esti_D_NEON$Estimator)
plot(esti_D_NEON$Observed, esti_D_NEON$Estimator)

##### Extract sample coverage ##### 
sampCover_NEON <- getSampleCoverage(iNEXT_object = sampleCompletenessNEON$objects, 
                                   sites = sites)

sampCover_NEON %>% 
  #filter(order == "1") %>% 
  filter(SC < 0.9) %>% 
  view()
