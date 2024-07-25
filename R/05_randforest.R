library(tidyverse)
library(ranger)
library(vip)


dat <- readRDS("./output/data_frame_for_modeling.RDS")
names(dat)
predictors <- c("total_developed_km2","total_impervious_km2","unclassified","primary_road","tertiary_road",
                "non_road_non_energy_impervious","secondary_road","lcmap_impervious","area_km2","open_water",
                "developed_open_space","developed_low_intensity","developed_medium_intensity",
                "developed_high_intensity","barren_land","deciduous_forest","evergreen_forest",
                "mixed_forest","shrub_scrub","herbaceous","hay_pasture","woody_wetlands",
                "emergent_herbaceous_wetlands","cultivated_crops","proportion_developed",
                "proportion_impervious","proportion_developed_high","proportion_developed_med",
                "proportion_developed_low","proportion_developed_openspace","proportion_cultivated_crops")
scaled_predictors <- c("proportion_developed",
                       "proportion_impervious","proportion_developed_high","proportion_developed_med",
                       "proportion_developed_low","proportion_developed_openspace","proportion_cultivated_crops")
outcomes <- c("proportion_unique_taxa","proportion_rare_taxa","proportion_pathogen",
              "proportion_mutualist","proportion_saprotroph","Observed","Shannon")


mod_list <- list()
for(i in seq_along(outcomes)){
  # build df with only one outcome
  new_df <- dat %>% dplyr::select(outcomes[i],all_of(scaled_predictors))
  # run full model on that
  full_mod <- glm(data=new_df,
                  formula = as.formula(paste(outcomes[i],"~",paste(scaled_predictors,collapse = " + "))))
  # reduce model
  step <- MASS::stepAIC(full_mod)
  step$formula
  # re-run model
  reduced_mod <- glm(data=new_df,
                     formula = step$formula)
  
  # save model
  current_mod <- broom::tidy(reduced_mod) %>% 
    mutate(outcome = outcomes[i],
           Rsquared = rsq::rsq(reduced_mod))
  # add to list
  mod_list[[i]] <- current_mod
  
}
mod_list %>% 
  purrr::reduce(full_join) %>% 
  mutate(across(is.numeric,function(x){round(x,4)})) %>% 
  write_csv("./output/reduced_model_output.csv")


# rare taxa drive richness/shannon?
dat %>% 
  ggplot(aes(x=proportion_rare_taxa,y=Observed)) +
  geom_point()
cor.test(dat$proportion_rare_taxa,dat$Observed)
cor.test(dat$proportion_rare_taxa,dat$Shannon)

ps <- readRDS("./output/ps_cleaned.RDS")
library(phyloseq)
core_asvs <- read_csv("./output/core_taxa_relabund.csv")[1:10,] %>% pluck("ASV")
core_names <- c("Plectosphaerella cucumerina",
                "Lemonniera centrosphaera",
                "Neopyrenochaeta annellidica",
                "Neonectria candida",
                "Nectriaceae sp.",
                "Cladosporium subuliforme",
                "Pleosporales sp.",
                "Neoascochyta desmazieri",
                "Pseudohalonectria adversaria",
                "Malasseziaceae sp.")


core_ps <- ps %>% 
  merge_samples(group = "location2") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  subset_taxa(taxa_names(ps) %in% core_asvs)
core_ps@tax_table[,7] <- core_names

melt <- psmelt(core_ps)
names(melt)
melt$Sample

new.pal <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")
melt %>% 
  mutate(Species = factor(Species, levels = core_names)) %>% 
  ggplot(aes(x=factor(Sample,levels=as.character(1:10)),
             y=Abundance,
             fill=Species)) +
  geom_col() +
  theme_bw() +
  labs(x="Location",
       y="Relative abundance",
       fill="Core microbiome\ntaxa") +
  theme(legend.title = element_text(face='bold',size = 18),
        legend.text = element_text(face='bold.italic',size=12),
        axis.title = element_text(face='bold',size=18),
        axis.text = element_text(face='bold',size=14)) +
  scale_fill_manual(values=new.pal)
ggsave("./output/figs/manuscript_figs/Figure_2_core_taxa.tiff",height = 5,width = 8,dpi=400)
ggsave("./output/figs/manuscript_figs/Figure_2_core_taxa.png",height = 5,width = 8,dpi=400)

ps %>% 
  transform_sample_counts(function(x){x/sum(x)})
melt2 <- ps %>% 
  merge_samples(group = "location2") %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt()

s <- 
melt2 %>% group_by(OTU) %>% 
  summarize(min = min(Abundance),
            max = max(Abundance)) %>% 
  arrange(desc(max))




melt2 %>% 
  # dplyr::filter(OTU %in% 
  #                 (melt2 %>% group_by(OTU) %>% 
  #                 summarize(mean = mean(Abundance)) %>% 
  #                 arrange(desc(mean)) %>% 
  #                 dplyr::filter(mean > 0.001) %>% 
  #                 pluck("OTU"))) %>% 
  ggplot(aes(x=Abundance)) +
  theme_minimal() + theme(strip.text = element_blank(),
                          strip.background = element_blank(),
                          axis.title = element_blank(),
                          axis.text = element_blank()) +
  geom_density(color='black') +
  facet_wrap(~OTU,scales = 'free') 
ggsave("./output/figs/manuscript_figs/SI_Fig_1_relabund_distributions.tiff",width = 12,height = 8,dpi=400)

s %>% 
  mutate(OTU = factor(OTU,levels=s$OTU)) %>% 
  ggplot(aes(x=OTU,ymin=min,ymax=max)) +
  geom_errorbar() +
  theme_void()


melt3 <- 
ps %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt()
rare_taxa_names <- taxa_names(ps)[which(taxa_sums(ps) / sum(sample_sums(ps)) < 0.001)]

n_rare_df <- 
melt3 %>% 
  mutate(rare = case_when(OTU %in% rare_taxa_names & Abundance > 0 ~ TRUE,
                          TRUE ~ FALSE)) %>% 
  group_by(Sample) %>% 
  summarize(N_rare = sum(rare)) %>% 
  mutate(location = str_remove(Sample,"A|B|C") %>% str_split("_") %>% map_chr(2) %>% factor(levels = as.character(1:10)))

n_rare_df %>% 
  ggplot(aes(x=location,y=N_rare)) +
  geom_boxplot() +
  labs(y="N rare taxa",x="Location") +
  theme_minimal() +
  theme(legend.title = element_text(face='bold',size = 18),
        legend.text = element_text(face='bold.italic',size=12),
        axis.title = element_text(face='bold',size=18),
        axis.text = element_text(face='bold',size=14))
ggsave("./output/figs/manuscript_figs/SI_Fig_2_N_rare_taxa_by_Location.tiff",height = 6,width = 6,dpi=400)

n_rare_df %>% 
  glm(data=.,
      formula = N_rare ~ as.numeric(location)) %>% 
  broom::tidy()

n_rare_df %>% arrange(desc(N_rare)) %>% tail

# rare taxa
rare_mod <- 
dat %>% 
  dplyr::select(proportion_pathogen,all_of(scaled_predictors)) %>% 
  ranger::ranger(data = .,
                 formula = proportion_pathogen ~ .,importance = 'permutation')
vip(rare_mod)
preds <- predict(rare_mod,dat)
plot(dat$proportion_cultivated_crops,preds$predictions)




