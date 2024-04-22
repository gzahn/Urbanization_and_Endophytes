# Make figures and tables for manuscript ####

## Setup ####
library(tidyverse)
library(patchwork)
library(broom)
library(VennDiagram)
library(phyloseq)
'%ni%' <- Negate("%in%")
### Load models ####
model_paths <- list.files("./output/models",full.names = TRUE, pattern = ".RDS")
for(i in model_paths){
  mod_name <- i %>% basename() %>% str_remove(".RDS")
  assign(mod_name,readRDS(i),envir = .GlobalEnv)
}
for(i in ls(pattern = "_mod$")){
  x <- get(i)
  print(formula(x))
}
sink(NULL)
sink("./output/models/model_specifications.txt")
print("gamma_diversity");gamma_div_mod %>% formula;gamma_div_mod %>% summary
print("mutualist_mod");mutualist_mod %>% formula;mutualist_mod %>% summary
print("pathogen_mod");pathogen_mod %>% formula;pathogen_mod %>% summary
print("permanova_mod");permanova_mod
print("raretaxa_mod");raretaxa_mod %>% formula;raretaxa_mod %>% summary
print("shannon_div_mod");shannon_div_mod %>% formula;shannon_div_mod %>% summary
print("uniquetaxa_mod");uniquetaxa_mod %>% formula;uniquetaxa_mod %>% summary
sink(NULL)


# core taxa and overlap between urban and rural
ps <- readRDS("./output/ps_cleaned.RDS")
rural <- ps %>% subset_samples(location < 5)
urban <- ps %>% subset_samples(location > 4)
rural_core <- microbiome::core_members(rural,detection = .1, prevalence = .5) %>% corncob::otu_to_taxonomy(data=rural)
urban_core <- microbiome::core_members(urban,detection = .1, prevalence = .5) %>% corncob::otu_to_taxonomy(data=urban)

rural <- rural %>% subset_taxa(taxa_sums(rural) > 0)
urban <- urban %>% subset_taxa(taxa_sums(urban) > 0)

set1 <- taxa_names(rural)
set2 <- taxa_names(urban)

venn.diagram(
  x = list(set1,set2),
  category.names = c("Rural sites","Urban sites"),
  filename = "./output/figs/venn_diagram.tiff",
  output = TRUE,
  height = 6,
  width = 6,units = "in",resolution = 400,imagetype = 'tiff',disable.logging = TRUE,        
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cex=1.5,fontface='bold',
  fill = c("darkblue","darkorange")
)

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

kableExtra::save_kable(kb,file = "./output/models/models_outputs_kable_table.png")




## Regression plots ####
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
        axis.text = element_text(face='bold',size=10),
        strip.background = element_rect(fill='white')) +
  labs(x="Proportion of impermeable surface in watershed",y="Response values")
ggsave("./output/figs/impermeable_surface_regression_plots.png",width = 6, height = 7,dpi=400)

(ps@tax_table[,2] %>% table) / ntaxa(ps)


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

