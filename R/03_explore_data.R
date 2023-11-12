
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(broom)
library(readxl)
library(corncob)
library(microbiome)
library(ecodist)
library(ggmap); packageVersion("ggmap")
source("./R/googlemap_styling.R")

set.seed(666)

# Load google maps API key from .Renviron 
ggmap::register_google(key = Sys.getenv("APIKEY")) # Key kept private

readRenviron("~/.Renviron")
options(warn = 1)
source("./R/plot_bar2.R")
set.seed(666)

############

ps <- readRDS("./output/ps_not-cleaned_no-metadata.RDS")
# add lat/lon metadata
meta <- read_xlsx("./data/aquatic-plant-proj-metadata.xlsx") %>% 
  mutate(sampleid = paste0("JM_",site_ID)) %>% 
  dplyr::filter(sampleid %in% ps@sam_data$sampleid) %>% 
  mutate(lat=lat %>% str_remove("_N") %>% as.numeric(),
         lon=lon %>% str_remove("_W") %>% as.numeric() * -1)
glimpse(meta)
meta

sam <- 
sample_data(
  full_join(
    microbiome::meta(ps), 
    meta
  )
)
sample_names(sam) <- meta$sampleid
ps <- phyloseq(sam,otu_table(ps),tax_table(ps))
# clean up taxa
# keep only unambiguous fungi
ps <- 
  ps %>% 
  subset_taxa(Kingdom == "k__Fungi") %>% 
  subset_taxa(!is.na(Phylum))

# convert location to factor and numeric 
ps@sam_data$location2 <- factor(ps@sam_data$location,
                                levels = as.character(c(1,2,3,4,5,6,9,10,7,8))) # actual river position order (downstream)
ps@sam_data$location <- ps@sam_data$location %>% as.numeric()


# export table of taxonomy and total read counts
data.frame(taxonomy=unname(corncob::otu_to_taxonomy(taxa_names(ps),ps)),
           total_reads=unname(taxa_sums(ps))) %>% 
  write_csv("./output/taxonomy_counts.csv")
  
plot_bar2(ps %>% transform_sample_counts(function(x){x/sum(x)}),fill = "Phylum")


# alpha diversity
diversity <- estimate_richness(ps,measures = c("Observed","Shannon")) %>% 
  mutate(location = ps@sam_data$location)
plot_richness(ps,measures = c("Observed"),x = "location2") + 
  geom_boxplot(fill="#4d6059") +
  labs(y="Fungal richness",x="Downstream position") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  scale_x_discrete(labels=as.character(1:10)) +
  geom_vline(xintercept = 5.5,linetype=2,color='red') +
  annotate('text',label="Edge of 'urban zone'",x=7,y=75)
ggsave("./output/figs/ASV_Richness.png",height = 4, width = 8,dpi=300)

ps@sam_data %>% 
  data.frame() %>% 
  ggplot(aes(x=lon,y=lat)) +
  geom_point()

# model alpha diversity
richness_mod <- diversity %>% 
  mutate(location=as.character(location) %>% as.numeric) %>% 
  glm(data=., formula = Observed ~ location)
summary(richness_mod)
tidy(richness_mod) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic()


# which taxa are dropping out as urbanization increases?
# do this in relative abundance
ps_ra <- 
ps %>% 
  transform_sample_counts(function(x){x/sum(x)})
sample_sums(ps_ra)

diversity %>% 
  mutate(location2=as.character(location) %>% as.numeric()) %>% 
  ggplot(aes(x=Observed,y=Shannon)) +
  geom_point() +
  geom_smooth(method='lm')

plot_richness(ps_ra,measures = c("Shannon"),x = "location2") + 
  geom_boxplot() +
  labs(y="Species Richness",x="Location") +
  theme_bw() +
  theme(strip.text = element_blank())


diversity %>%
  mutate(lat=ps_ra@sam_data$lat,lon=ps_ra@sam_data$lon) %>% 
  mutate(location2=as.character(location) %>% as.numeric()) %>% 
  ggplot(aes(x=lat,y=lon,color=Observed)) +
  geom_point(alpha=.5,size=4) +
  scale_color_viridis_c() +
  theme_bw()


# Multiple regression on matrices
gps_dist <- vegan::vegdist(microbiome::meta(ps_ra) %>% select(lat,lon),
                    method = "euclidean")
asv_dist <- vegan::vegdist(otu_table(ps_ra),method = "bray") 

mrm <- MRM(asv_dist ~ gps_dist)
data.frame(mrm) %>% 
  write_csv("./output/MRM_stats_table.csv")

plot_bar2(ps_ra,fill="Phylum") +
  scale_fill_viridis_d()
ggsave("./output/figs/stacked_barchart_Phylum.png",dpi=300,height = 4,width = 4)

# Map of sites ####
meta <- microbiome::meta(ps)

mapstyle <- rjson::fromJSON(file = "./R/mapstyle3.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

area <- 
  ggmap::get_googlemap(center = c(lon = ps@sam_data$lon %>% mean, 
                                  lat = ps@sam_data$lat %>% mean),
                       zoom = 12,
                       scale = 2,
                       style=mapstyle)

ggmap::ggmap(area) +
  geom_point(data=meta,aes(x=lon,y=lat),color='darkorange',size=4,alpha=.75) 
ggsave("./output/figs/map.png",height = 8,width = 8,dpi=300)

# try looking for taxa tables (existing taxa only) for each location
# make list of all taxa names present in each location
taxa_list <- list()
# i <- 1
for(i in 1:10){
  ps_x <- 
    ps %>% 
    subset_samples(location == i)
  print(ps_x)
  taxa_x <- 
    ps_x %>% 
    # merge_samples("location",fun = "sum") %>% 
    subset_taxa(taxa_sums(ps_x) > 0) %>% 
    taxa_names() %>% 
    as.data.frame()
  names(taxa_x) <- as.character(i)
  df <- taxa_x %>% table %>% data.frame()
  names(df) <- c("ASV","Freq")
  df$Location <- as.character(i)
  taxa_list[[i]] <- df
}
# build data frame from results
names(taxa_list) <- as.character(1:10)
df <- taxa_list[] %>% 
  purrr::reduce(full_join)

# subset to non-duplicated taxa
df[!(df$ASV %>% duplicated()),] %>% 
  group_by(Location) %>% 
  summarize(rare_count=sum(Freq)) %>% 
  mutate(Location=as.numeric(Location) %>% factor(levels=as.character(c(1,2,3,4,5,6,9,10,7,8)))) %>% 
  arrange(Location) %>% 
  ggplot(aes(x=Location,y=rare_count)) +
  geom_col(fill="#4d6059") +
  labs(y="N unique taxa",
       title = "Number of unique taxa in each location",
       x="Downstream position") +
  scale_x_discrete(labels=as.character(1:10)) +
  geom_vline(xintercept = 5.5,linetype=2,color='red') +
  annotate('text',label="Edge of 'urban zone'",x=7,y=75) +
  theme_bw()
ggsave("./output/figs/N_unique_taxa_at_each_location.png", width = 8, height = 4,dpi=300)



# taxonomy of taxa unique to only one location
unique_asvs <- df[!(df$ASV %>% duplicated()),"ASV"] %>% as.character()

unique_asvs_taxonomy <- 
ps %>% 
  subset_taxa(taxa_names(ps) %in% unique_asvs) %>% 
  tax_table()
x <- 
data.frame(unique_asvs_taxonomy) %>% 
  dplyr::filter(!is.na(Genus)) %>% 
  pluck("Genus") %>% 
  str_remove("g__") %>% 
  table() %>% 
  data.frame() %>% 
  arrange((Freq))
names(x) <- c("Genus","Freq")
x %>% 
  dplyr::filter(Freq > 1) %>% 
  mutate(Genus=factor(Genus,levels=x$Genus)) %>% 
  ggplot(aes(y=Genus,x=Freq)) +
  geom_col() +
  theme_bw()


ps_comp <- microbiome::transform(ps, "compositional")

# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(ps_comp, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

core_taxa <- core_members(ps_comp, detection = 0, prevalence = 2/3) %>% otu_to_taxonomy(data=ps_comp)

data.frame(CoreTaxa=core_taxa %>% 
             str_remove_all(".__") %>% 
             str_remove("Fungi_")) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic()



# quick ordination
set.seed(666)
ord <- ordinate(ps_comp,method = "NMDS")
plot_ordination(ps_comp,ord,color = "location") + 
  scale_color_viridis_c(end=.9,option = "B") +
  theme_minimal()







# quick and dirty diffabund analysis

library(corncob)
source("./R/bbdml_helper.R")
ps@sam_data %>% str
da_analysis <- differentialTest(formula = ~ location2, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", boot = FALSE,
                                         data = ps,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)
plot(da_analysis) +
  theme(legend.position = 'none')

taxa_names(ps)
da_analysis$significant_taxa %>% otu_to_taxonomy(data=ps)
bbd <- bbdml(formula = CTGAGTACTACACTCTCTACCCTTTGTGAACTATTATACCTGTTGCTTCGGCGGCGCCCGCAAGGGTGCCCGCCGGTCTCATCAGAATCTCTGTTTTCGAACCCGACGATACTTCTGAGTGTTCTTAGCGAACTGTCA ~ location2,
      phi.formula = ~ location2,
      data = ps)

plot(bbd, color="location2")
plot_multi_bbdml


