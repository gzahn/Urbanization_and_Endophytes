# Make figures and tables for manuscript ####

## Setup ####
library(tidyverse)
library(patchwork)
library(broom)

### Load models ####
model_paths <- list.files("./output/models",full.names = TRUE, pattern = ".RDS")
for(i in model_paths){
  mod_name <- i %>% basename() %>% str_remove(".RDS")
  assign(mod_name,readRDS(i),envir = .GlobalEnv)
}

### Load metadata ####
meta <- readRDS("./output/cleaned_full_metadata.RDS")

### Load plots ####
diffabund_plot <- readRDS("./output/figs/diffabund_plot_high_development_categorical.RDS")
fig.paths <- list.files(pattern = "png",full.names = TRUE,recursive = TRUE)


### Combine models ####


richness_mod <- gamma_div_mod
raretaxa_mod
shannon_div_mod
options(scipen = 999)
kb <- broom::tidy(raretaxa_mod) %>% mutate(outcome="Rare taxa") %>% 
  full_join(broom::tidy(richness_mod) %>% mutate(outcome="Richness")) %>% 
  full_join(broom::tidy(shannon_div_mod) %>% mutate(outcome="Shannon div.")) %>% 
  mutate(across(where(is.numeric),function(x){round(x,2)})) %>% 
  mutate(p.value = kableExtra::cell_spec(p.value, bold = ifelse(p.value < 0.05,TRUE,FALSE))) %>% 
  kableExtra::kable(escape = FALSE) %>% 
  kableExtra::kable_classic() %>% 
  kableExtra::row_spec(0,bold=TRUE)
kb

# kableExtra::save_kable(kb,file = "./output/models/kable_table.png")




## Regression plots ####
meta %>% names

simple <- 
meta %>% 
  dplyr::select(Shannon,Observed,richness_proportion,proportion_unique_taxa,proportion_rare_taxa,
                proportion_cultivated_crops,proportion_developed,proportion_developed_high,proportion_impervious)

simple %>% 
  dplyr::select(proportion_cultivated_crops,proportion_developed,proportion_developed_high,proportion_impervious) %>% 
  GGally::ggpairs()

simple %>% 
  pivot_longer(all_of(c("Shannon","Observed","richness_proportion","proportion_unique_taxa","proportion_rare_taxa")),
               names_to = "outcome") %>% 
  mutate(outcome = case_when(outcome == "Shannon" ~ "Shannon div.",
                             outcome == "Observed" ~ "Richness",
                             outcome == "richness_proportion" ~ "Local gamma div.",
                             outcome == "proportion_unique_taxa" ~ "Unique taxa proportion",
                             outcome == "proportion_rare_taxa" ~ "Rare taxa proportion")) %>%
  dplyr::filter(outcome != "Unique taxa proportion" & outcome != "Local gamma div.") %>% 
  ggplot(aes(x=proportion_impervious,y=value)) +
  geom_point(size=3) +
  geom_smooth(method='lm',se=FALSE, color="darkgreen") +
  facet_wrap(~outcome,scales = 'free',ncol = 1) +
  theme_bw() +
  theme(strip.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=14),
        axis.text = element_text(face='bold',size=10)) +
  labs(x="Proportion of impermeable surface in watershed",y="Response values")
ggsave("./output/figs/impermeable_surface_regression_plots.png",width = 6, height = 7,dpi=400)







# examine all significant predictors
simple_long <- 
simple %>% 
  pivot_longer(all_of(c("proportion_cultivated_crops","proportion_developed","proportion_developed_high","proportion_impervious")),
               names_to = "predictor",values_to = "predictor_value")

simple_long %>% 
  ggplot(aes(x=predictor_value,y=Shannon)) +
  geom_point() +
  facet_wrap(~predictor,scales = 'free') +
  theme_minimal()

simple_long %>% 
  ggplot(aes(x=predictor_value,y=Observed)) +
  geom_point() +
  facet_wrap(~predictor,scales = 'free') +
  theme_minimal()

simple_long %>% 
  ggplot(aes(x=predictor_value,y=richness_proportion)) +
  geom_point() +
  facet_wrap(~predictor,scales = 'free') +
  theme_minimal()

simple_long %>% 
  ggplot(aes(x=predictor_value,y=proportion_unique_taxa)) +
  geom_point() +
  facet_wrap(~predictor,scales = 'free') +
  theme_minimal()

simple_long %>% 
  ggplot(aes(x=predictor_value,y=proportion_rare_taxa)) +
  geom_point() +
  facet_wrap(~predictor,scales = 'free') +
  theme_minimal()

