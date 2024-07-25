# Make figures and tables for manuscript ####

## Setup ####
library(tidyverse)
library(patchwork)
library(broom)
library(VennDiagram)
library(phyloseq)
library(zahntools)
'%ni%' <- Negate("%in%")
ra <- function(x){x/sum(x)}
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


# refit with location as covariate to make sure this isn't just spatial!
raretaxa_mod$formula
full <- glm(data = raretaxa_mod$data,
    formula = proportion_rare_taxa ~ proportion_developed + proportion_impervious + 
      proportion_developed_high + proportion_developed_med + proportion_developed_low + location) 
step <- MASS::stepAIC(full)
glm(data = raretaxa_mod$data,
    formula = step$formula) %>% 
  summary



gamma_div_mod$formula
glm(data = gamma_div_mod$data,
    formula = richness_proportion ~ proportion_impervious + proportion_developed_low + 
      proportion_cultivated_crops + location) %>% 
  summary
options(scipen = 999)

shannon_div_mod$formula
glm(data = shannon_div_mod$data,
    formula = Shannon ~ proportion_developed + proportion_developed_low + proportion_cultivated_crops + location) %>% 
  summary

mutualist_mod$formula
glm(data = mutualist_mod$data,
    formula = proportion_mutualist ~ proportion_developed + proportion_impervious + location) %>% 
  summary

glm(data = pathogen_mod$data,
    formula = proportion_pathogen ~ proportion_developed + proportion_developed_high + location) %>% 
  summary

# core taxa and overlap between urban and rural
ps <- readRDS("./output/ps_cleaned.RDS")
core <- microbiome::core_members(ps,detection = .1, prevalence = .5) %>% corncob::otu_to_taxonomy(data=ps)
rural <- ps %>% subset_samples(location < 5)
urban <- ps %>% subset_samples(location > 4)
rural_core <- microbiome::core_members(rural,detection = .1, prevalence = .5) %>% corncob::otu_to_taxonomy(data=rural)
urban_core <- microbiome::core_members(urban,detection = .1, prevalence = .5) %>% corncob::otu_to_taxonomy(data=urban)

ps_core <- ps %>% 
  transform_sample_counts(ra) %>% 
  subset_taxa(taxa_names(ps) %in% names(core))
ps_core %>% otu_table() %>% as("matrix") %>% heatmap

((ps_core %>% taxa_sums()) / nsamples(ps_core)) %>% unname


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
  ggplot(aes(x=proportion_developed,y=value)) +
  geom_point(size=3) +
  geom_smooth(method='lm',se=FALSE, color="darkgreen") +
  facet_wrap(~outcome,scales = 'free',ncol = 1) +
  theme_bw() +
  theme(strip.text = element_text(face='bold',size=18),
        axis.title = element_text(face='bold',size=18),
        axis.text = element_text(face='bold',size=14),
        strip.background = element_blank()) +
  labs(x="Proportion of developed land\nsurface in watershed",y="Response values")
ggsave("./output/figs/impermeable_surface_regression_plots.png",width = 6, height = 8,dpi=400)
ggsave("./output/figs/impermeable_surface_regression_plots.pdf",width = 6, height = 8,dpi=600)
ggsave("./output/figs/manuscript_figs/Figure_3_regressions.pdf",width = 6, height = 8,dpi=600)


ps@sam_data$location2 <- factor(ps@sam_data$location,levels = as.character(1:10))
ps

  merge_samples(group = 'location2') %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(x = location,fill="Class")

(ps@tax_table[,2] %>% table) / ntaxa(ps)

citation("MASS")
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

# guilds vs urbanization plots
df <- shannon_div_mod$data

df %>% 
  mutate(impervious = total_impervious_km2 / area_km2,
         total_road_perc = sum(primary_road, secondary_road,tertiary_road) / area_km2,
         ) %>% 
  pivot_longer(c(proportion_pathogen,proportion_mutualist),
               names_to = "guild",
               names_prefix = "proportion_",
               values_to = "Proportion") %>% 
ggplot(aes(x=location,y=Proportion,color=guild)) +
  geom_point() +
  geom_smooth(method='lm')

df %>% 
  ggplot(aes(x=factor(location),y=proportion_pathogen)) +
  geom_boxplot() +
  geom_smooth()

df %>% names
# land use plots
imperm_vars <- c("primary_road","secondary_road","tertiary_road","non_road_non_energy_impervious",
                    "lcmap_impervious","total_impervious_km2")
developed_vars <- c("developed_low_intensity", "developed_medium_intensity",
                    "developed_high_intensity","developed_open_space","total_developed_km2")
permeable_vars <- c("open_water","barren_land","deciduous_forest","evergreen_forest","mixed_forest",
                      "shrub_scrub","herbaceous","hay_pasture","woody_wetlands",
                      "emergent_herbaceous_wetlands","cultivated_crops")

long <- 
df %>% 
  pivot_longer(all_of(imperm_vars),
               names_to = "impermeable_surface_category",
               values_to = "impermeable_surface_area_km2") %>% 
  pivot_longer(all_of(developed_vars),
               names_to = "development_category",
               values_to = "developed_area_km2") %>% 
  pivot_longer(all_of(permeable_vars),
               names_to = "permeable_surface_category",
               values_to = "permeable_surface_area_km2") %>% 
  mutate(impermeable_surface_category = impermeable_surface_category %>% str_replace_all("_"," ") %>% str_remove(" km2"),
         development_category = development_category %>% str_replace_all("_"," ") %>% str_remove(" km2"),
         permeable_surface_category = permeable_surface_category %>% str_replace_all("_"," ") %>% str_remove(" km2"))

# long %>% 
#   mutate(location = factor(location,levels=c(as.character(1:10))))


# permeable surface plot
permeable_surface_plot <- 
long %>% 
  ggplot(aes(x=location,y=permeable_surface_area_km2,color=permeable_surface_category)) +
  geom_point(size=3) +
  geom_path() +
  scale_color_viridis_d(end=.85) +
  theme_bw() +
  labs(x="Sample location",
       y="Surface area (km2)",
       color="Permeable surface\ncategory") +
  theme(axis.title = element_text(face='bold',size=18),
        axis.text = element_text(face='bold',size=14),
        legend.title = element_text(face='bold',size=18),
        legend.text = element_text(face='bold',size=12),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks=1:10)

# impermeable surface plot
impermeable_surface_plot <- 
  long %>% 
  ggplot(aes(x=location,y=impermeable_surface_area_km2,color=impermeable_surface_category)) +
  geom_point(size=3) +
  geom_path() +
  scale_color_viridis_d(end=.85) +
  theme_bw() +
  labs(x="Sample location",
       y="Surface area (km2)",
       color="Impermeable surface\ncategory") +
  theme(axis.title = element_text(face='bold',size=18),
        axis.text = element_text(face='bold',size=14),
        legend.title = element_text(face='bold',size=18),
        legend.text = element_text(face='bold',size=12)) +
  scale_x_continuous(breaks=1:10)

# developed surface plot
developed_surface_plot <- 
  long %>% 
  ggplot(aes(x=location,y=developed_area_km2,color=development_category)) +
  geom_point(size=3) +
  geom_path() +
  scale_color_viridis_d(end=.85) +
  theme_bw() +
  labs(x="Sample location",
       y="Surface area (km2)",
       color="Development category") +
  theme(axis.title = element_text(face='bold',size=18),
        axis.text = element_text(face='bold',size=14),
        legend.title = element_text(face='bold',size=18),
        legend.text = element_text(face='bold',size=12),
        axis.title.x = element_blank()) +
  scale_x_continuous(breaks=1:10)
developed_surface_plot



# patch together into one plot
patch <-
permeable_surface_plot / developed_surface_plot / impermeable_surface_plot & ylab(NULL)

patch2 <- 
  wrap_elements(panel = patch) +
  labs(tag = "Total area (km2)") +
  theme(
    plot.tag = element_text(size = 18,face='bold',angle=90),
    plot.tag.position = "left"
  )
patch2 
ggsave(plot = patch2,filename = "./output/figs/manuscript_figs/land_use.pdf",height = 10,width = 8,dpi=600,device = 'pdf')
ggsave(plot = patch2,filename = "./output/figs/manuscript_figs/land_use.png",height = 10,width = 8,dpi=600,device = 'png')
