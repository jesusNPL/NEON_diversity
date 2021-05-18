
demon_TRAIT_summarize <- function(df, trait_name_IN, trait_name_TO, tax_level) {
  
  tt <- df %>% 
    filter(trait_name == trait_name_IN) 
  
  tt <- transform(tt, trait_value = as.numeric(trait_value))
  
  if (tax_level == "species") {
    t_summ <- tt %>% 
      group_by(SciNames) %>% 
      summarise(mean_val = mean(trait_value, na.rm = TRUE),  
                sd_val = sd(trait_value)) 
    t_summ$Trait <- trait_name_TO 
    t_summ$Tax_level <- tax_level
  } else if (tax_level == "genus") {
    t_summ <- tt %>% 
      group_by(Genus) %>% 
      summarise(mean_val = mean(trait_value, na.rm = TRUE),  
                sd_val = sd(trait_value)) 
    t_summ$Trait <- trait_name_TO 
    t_summ$Tax_level <- tax_level 
  } else if (tax_level == "family") {
    t_summ <- tt %>% 
      group_by(Family) %>% 
      summarise(mean_val = mean(trait_value, na.rm = TRUE),  
                sd_val = sd(trait_value)) 
    t_summ$Trait <- trait_name_TO 
    t_summ$Tax_level <- tax_level
  }
  
  names(t_summ) <- c(names(t_summ)[1], paste0("mean_", trait_name_TO), paste0("SD_", trait_name_TO), "Trait", "Tax_level")
  return(t_summ)
}

# xx <- demon_TRAIT_summarize(df = traits_NEON, trait_name_IN = "whole plant height", trait_name_TO = "WPH", tax_level = "family")

# This function summarize the trait data from BIEN genus level only 

demon_TRAITgenus_summarize <- function(df, trait_name_IN, trait_name_TO) {
  tt <- df %>% 
    filter(trait_name == trait_name_IN) 
  
  tt <- transform(tt, trait_value = as.numeric(trait_value))
  
  if (tax_level == "species") {
    t_summ <- tt %>% 
      group_by(SciNames) %>% 
      summarise(mean_val = mean(trait_value, na.rm = TRUE),  
                sd_val = sd(trait_value)) 
    t_summ$Trait <- trait_name_TO 
    t_summ$Tax_level <- tax_level
  } else if (tax_level == "genus") {
    t_summ <- tt %>% 
      group_by(Genus) %>% 
      summarise(mean_val = mean(trait_value, na.rm = TRUE),  
                sd_val = sd(trait_value)) 
    t_summ$Trait <- trait_name_TO 
    t_summ$Tax_level <- tax_level 
  } else if (tax_level == "family") {
    t_summ <- tt %>% 
      group_by(Family) %>% 
      summarise(mean_val = mean(trait_value, na.rm = TRUE),  
                sd_val = sd(trait_value)) 
    t_summ$Trait <- trait_name_TO 
    t_summ$Tax_level <- tax_level
  }
  
  names(t_summ) <- c(names(t_summ)[1], paste0("mean_", trait_name_TO), paste0("SD_", trait_name_TO), "Trait", "Tax_level")
  return(t_summ)
  
}