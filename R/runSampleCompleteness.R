library(tidyverse)
library(ggdist)
theme_set(bayesplot::theme_default(base_family = "sans"))

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
  filter(Diversity == "Species richness") %>% 
  mutate(Difference1 = (Estimator - Observed), Difference2 = (Observed - Estimator))

esti_SR_NEON %>% 
  filter(Difference1 < 100) %>% 
  ggplot(aes(x = Site, y = Difference1, colour = Site)) + 
  #geom_boxplot(outlier.shape = NA, coef = 0, width = 0.15) + 
  stat_gradientinterval(position = "dodge", #p_limits = c(0.025, 0.975), 
                        thickness = stat(0.5), 
                        justification = stat(0.25), fill = "darkgray", 
                        colour = "darkgray") + 
  coord_flip()
  

view(esti_SR_NEON)

cor.test(esti_SR_NEON$Observed, esti_SR_NEON$Estimator)

esti_SR_NEON %>% 
  filter(Estimator < 80) %>% 
  ggplot(aes(x = Estimator, y = Observed)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_cartesian(xlim = c(0, 70), ylim = c(0, 70))

obs <- esti_SR_NEON %>% 
  group_by(Site) %>% 
  filter(Site != "JERC") %>% 
  summarise(Observed = mean(Observed), Rarefied = mean(Estimator)) %>% 
  select(Site, Observed) %>% 
  mutate(SR = "Observed", paired = 2:length(unique(esti_SR_NEON$Site))) %>% 
  rename(value = Observed)

rare <- esti_SR_NEON %>% 
  group_by(Site) %>% 
  filter(Site != "JERC") %>% 
  summarise(Observed = mean(Observed), Rarefied = mean(Estimator)) %>% 
  select(Site, Rarefied) %>% 
  mutate(SR = "Estimated", paired = 2:length(unique(esti_SR_NEON$Site))) %>% 
  rename(value = Rarefied)

dat_SR <- rbind(obs, rare) 

dat_SR %>%  
  #ggplot(aes(x = value, y = reorder(Site, value))) +
  ggplot(aes(x = value, y = Site)) + 
  geom_line(aes(group = paired), color = "black", alpha = 0.7, size = 1.5)+
  geom_point(aes(color = SR), alpha = 0.7, size = 7) +
  theme_classic(24) +
  theme(legend.position = "bottom") + 
  #scale_color_brewer(palette ="Accent", direction = -1) +
  scale_color_manual(values = c("red", "darkgreen")) + 
  labs(x = NULL, y = NULL)

## Shannon
esti_H_NEON <- sampleCompletenessNEON$Estimators %>% 
  filter(Diversity == "Shannon diversity")

cor.test(esti_H_NEON$Observed, esti_H_NEON$Estimator)
plot(esti_H_NEON$Observed, esti_H_NEON$Estimator)

esti_H_NEON %>% 
  ggplot(aes(x = Estimator, y = Observed)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 30))

obs <- esti_H_NEON %>% 
  group_by(Site) %>% 
  filter(Site != "JERC") %>% 
  summarise(Observed = mean(Observed), Rarefied = mean(Estimator)) %>% 
  select(Site, Observed) %>% 
  mutate(H = "Observed", paired = 2:length(unique(esti_H_NEON$Site))) %>% 
  rename(value = Observed)

rare <- esti_H_NEON %>% 
  group_by(Site) %>% 
  filter(Site != "JERC") %>% 
  summarise(Observed = mean(Observed), Rarefied = mean(Estimator)) %>% 
  select(Site, Rarefied) %>% 
  mutate(H = "Estimated", paired = 2:length(unique(esti_H_NEON$Site))) %>% 
  rename(value = Rarefied)

dat_H <- rbind(obs, rare) 

dat_H %>%  
  #ggplot(aes(x = value, y = reorder(Site, value))) +
  ggplot(aes(x = value, y = Site)) + 
  geom_line(aes(group = paired), color = "black", alpha = 0.7, size = 1.5)+
  geom_point(aes(color = H), alpha = 0.7, size = 7) +
  theme_classic(24) +
  theme(legend.position = "top") + 
  #scale_color_brewer(palette ="Accent", direction = -1) +
  scale_color_manual(values = c("red", "darkgreen")) + 
  labs(x = NULL, y = NULL)

## Simpson
esti_D_NEON <- sampleCompletenessNEON$Estimators %>% 
  filter(Diversity == "Simpson diversity")

cor.test(esti_D_NEON$Observed, esti_D_NEON$Estimator)
plot(esti_D_NEON$Observed, esti_D_NEON$Estimator)

esti_D_NEON %>% 
  ggplot(aes(x = Estimator, y = Observed)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 30))

obs <- esti_D_NEON %>% 
  group_by(Site) %>% 
  filter(Site != "JERC") %>% 
  summarise(Observed = mean(Observed), Rarefied = mean(Estimator)) %>% 
  select(Site, Observed) %>% 
  mutate(D = "Observed", paired = 2:length(unique(esti_H_NEON$Site))) %>% 
  rename(value = Observed)

rare <- esti_D_NEON %>% 
  group_by(Site) %>% 
  filter(Site != "JERC") %>% 
  summarise(Observed = mean(Observed), Rarefied = mean(Estimator)) %>% 
  select(Site, Rarefied) %>% 
  mutate(D = "Estimated", paired = 2:length(unique(esti_H_NEON$Site))) %>% 
  rename(value = Rarefied)

dat_D <- rbind(obs, rare) 

dat_D %>%  
  #ggplot(aes(x = value, y = reorder(Site, value))) +
  ggplot(aes(x = value, y = Site)) + 
  geom_line(aes(group = paired), color = "black", alpha = 0.7, size = 1.5)+
  geom_point(aes(color = D), alpha = 0.7, size = 7) +
  theme_classic(24) +
  theme(legend.position = "top") + 
  #scale_color_brewer(palette ="Accent", direction = -1) +
  scale_color_manual(values = c("red", "darkgreen")) + 
  labs(x = NULL, y = NULL)

##### Extract sample coverage #####  
load("output/NEON_sample_completeness.RData")

sites <- unique(sampleCompletenessNEON$Estimators$Site)

sampCover_NEON <- getSampleCoverage(iNEXT_object = sampleCompletenessNEON$objects, 
                                   sites = sites)

sampCover_NEON <- sampCover_NEON %>% 
  filter(Site != "JERC")

sampCover_NEON %>% 
  filter(order == "0") %>% 
  filter(SC < 0.9) %>% 
  count()

sampCover_NEON %>% 
  filter(order == "1") %>% 
  filter(SC > 0.9) %>% 
  count()

sampCover_NEON %>% 
  filter(order == "2") %>% 
  filter(SC <= 0.9) %>% 
  count()

round((923*100)/(nrow(sampCover_NEON)/3), 1)

barCOLS <- c("q2" = scales::alpha("red", 0.5), 
             "q1" = scales::alpha("darkblue", 0.5), 
             "q0" = scales::alpha("darkgray", 0.5))

sampCover_NEON$order <- factor(sampCover_NEON$order)

sampCover_NEON <- sampCover_NEON %>% 
  mutate(q = recode(order, "0" = "q0", "1" = "q1", "2" = "q2"))

pdf(file = "output/sampleCompleteness.pdf", 
    width = 5, height = 4) 

sampCover_NEON %>% 
  #filter(metric == "PCD") %>% 
  drop_na() %>% 
  ggplot(aes(x = q, 
             y = SC, 
             color = q, 
             fill = q)
  ) + 
  geom_jitter(shape = 16, position = position_jitter(0.15), alpha = 0.25) + 
  geom_boxplot(width = 0.3, notch = TRUE, notchwidth = 0.5) + 
  scale_color_manual(values = barCOLS, 
                     name = NULL
  ) + 
  geom_hline(yintercept = 0.9, col = "black", alpha = 0.50, lty = 2, size = 1.5) + 
  scale_fill_manual(values = barCOLS, 
                    name = NULL
  ) + 
  scale_y_continuous(n.breaks = 5, labels = scales::percent) + #, limits = c(0.8, 1.4)) + 
  labs(x = NULL, y = "Sample completeness") + 
  theme(
    legend.position = "none", 
    legend.title = element_text(size = 15), 
    legend.text = element_text(size = 15), 
    axis.text = element_text(size = 15, colour = "black"), 
    axis.text.x = element_text(hjust = 0.7), 
    axis.title = element_text(size = 20), 
    axis.line.y = element_line(size = 0.5), 
    axis.line.x = element_line(size = 0.5)
  )

dev.off()
