
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(broom)
library(readxl)
library(corncob)
library(microbiome)
library(ecodist)
library(ggmap); packageVersion("ggmap")
library(FUNGuildR)
library(fungaltraits)
library(DECIPHER)
library(phangorn)
library(Biostrings)

source("./R/bbdml_helper.R")
source("./R/googlemap_styling.R")
source("./R/palettes.R")

'%ni%' <- Negate('%in%')
set.seed(666)

# Load google maps API key from .Renviron 
ggmap::register_google(key = Sys.getenv("APIKEY")) # Key kept private

readRenviron("~/.Renviron")
options(warn = 1)
source("./R/plot_bar2.R")
set.seed(666)

# Load clean ps object ####

ps <- readRDS("./output/ps_not-cleaned_no-metadata.RDS")

# add lat/lon metadata
meta <- read_xlsx("./metadata/aquatic-plant-proj-metadata.xlsx") %>% 
  mutate(sampleid = paste0("JM_",site_ID)) %>% 
  dplyr::filter(sampleid %in% ps@sam_data$sampleid) %>% 
  mutate(lat=lat %>% str_remove("_N") %>% as.numeric(),
         lon=lon %>% str_remove("_W") %>% as.numeric())

sam <- 
sample_data(
  full_join(
    microbiome::meta(ps), 
    meta
  )
)
sample_names(sam) <- meta$sampleid
ps <- phyloseq(sam,otu_table(ps),tax_table(ps),phy_tree(ps))
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
  


# Alpha diversity ####
diversity <- estimate_richness(ps,measures = c("Observed","Shannon")) %>% 
  mutate(location = ps@sam_data$location)
plot_richness(ps,measures = c("Observed"),x = "location2") + 
  geom_boxplot(fill="#4d6059") +
  labs(y="Fungal richness",x="Downstream position") +
  theme_bw() +
  theme(strip.text = element_blank()) +
  scale_x_discrete(labels=as.character(1:10)) +
  geom_vline(xintercept = 5.5,linetype=2,color='red') +
  annotate('text',label="Edge of 'urban zone'",x=7,y=80) +
  geom_segment(aes(x=7,y=75,xend=5.5,yend=60),arrow = arrow(length = unit(.5,'cm')))
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

x <- ps_ra
x@tax_table[,2] <- x@tax_table[,2] %>% str_remove("p__")
meta <- microbiome::meta(ps)
sample_names(x) <- sample_names(x) %>% str_remove("JM_")
melt <- psmelt(x)
melt$sample_ID <- factor(str_remove(melt$sample_ID,"JM_"), levels = meta %>% 
                           arrange(desc(location)) %>% 
                           pluck('sample_ID') %>% 
                           str_remove("JM_"))
melt %>% 
  ggplot(aes(x=sample_ID,y=Abundance,fill=Phylum)) +
  geom_col(position = 'stack') +
# plot_bar2(x,fill="Phylum") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.title = element_text(size=14,face='bold'),
        legend.title = element_text(size=14,face='bold'),
        axis.text.x = element_text(angle=90,hjust=1,vjust=.5,face='bold')) +
  labs(x="Sample",y="Relative abundance")
ggsave("./output/figs/stacked_barchart_Phylum.png",dpi=300,height = 4,width = 6)

# Map of sites ####

mapstyle <- rjson::fromJSON(file = "./R/mapstyle3.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

area <- 
  ggmap::get_googlemap(center = c(lon = ps@sam_data$lon %>% mean, 
                                  lat = ps@sam_data$lat %>% mean),
                       zoom = 12,
                       scale = 2,
                       style=mapstyle)

ggmap::ggmap(area) +
  geom_point(data=meta,aes(x=lon,y=lat),color='darkorange',size=4) +
  geom_label(data=meta,aes(x=lon,y=lat,label=location),color='darkorange',size=4,alpha=.75)

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
  mutate(Location=as.numeric(Location) %>% factor(levels=as.character(1:10))) %>% 
  arrange(Location) %>% 
  ggplot(aes(x=Location,y=rare_count)) +
  geom_col(fill="#4d6059") +
  labs(y="N unique ASVs",
       # title = "Number of unique ASVs in each location",
       x="Downstream position") +
  scale_x_discrete(labels=as.character(1:10)) +
  geom_vline(xintercept = 5.5,linetype=2,color='red') +
  annotate('text',label="Edge of 'urban zone'",x=7,y=80) +
  geom_segment(aes(x=7,y=75,xend=5.5,yend=60),arrow = arrow(length = unit(.5,'cm'))) +
  theme_bw() +
  theme(axis.title = element_text(size=14,face='bold'),
        legend.title = element_text(size=14,face='bold'),
        axis.text.x = element_text(angle=90,hjust=1,vjust=.5,face='bold'))
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
  theme_bw() +
  theme(axis.title = element_text(size=14,face='bold'),
        legend.title = element_text(size=14,face='bold'),
        axis.text.x = element_text(size=14,face='bold')) +
  labs(x="Total observations of genus",subtitle = "Genera only observed at a single location")


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
  write_csv("./output/core_taxa_0_66.csv")

data.frame(CoreTaxa=core_taxa %>% 
             str_remove_all(".__") %>% 
             str_remove("Fungi_")) %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic()


# Beta diversity ####
# quick ordination
set.seed(666)
urban.sites <- c('10A','10B','10C','9A','9B','9C','8A','8B','8C',
                 '7A','7B','7C','6A','6B','6C','5A','5B','5C')
rural.sites <- c('4A','4B','4C','3A','3B','3C','2A','2B','2C','1A','1B','1C')

ps_comp@sam_data$urban <- ifelse(ps_comp@sam_data$site_ID %in% urban.sites,'urban','rural')
ord <- ordinate(ps_comp,method = "NMDS",distance = 'unifrac')
plot_ordination(ps_comp,ord,color = "urban") + 
  scale_color_viridis_d(end=.9,option = "B") +
  theme_minimal()
ps_comp@sam_data


# FunGuild analysis ####
db <- fungaltraits::fungal_traits()


guilds <- 
funguild_assign(data.frame(Taxonomy=paste(ps@tax_table[,1],
                                          ps@tax_table[,2],
                                          ps@tax_table[,3],
                                          ps@tax_table[,4],
                                          ps@tax_table[,5],
                                          ps@tax_table[,6],
                                          ps@tax_table[,7],
                                          sep=",")
))


guilds$guild

# add guild assignments to "pseudo-tax_table"
# (this replaces and renames the "Kingdom" column)
ps@tax_table[,1] <- guilds$guild
attributes(ps@tax_table)$dimnames[[2]][1] <- "Guild"
ps@tax_table[,1] %>% unique
# just using "mycorrhizal" as the keyword...
mutualist_guilds <- 
  grep("[M,m]ycorrhizal|[E,e]ndophyte",(ps@tax_table[,1]),value = TRUE) %>% 
  unique()

# identify "saprotrophs"
saprotroph_guilds <- 
  grep("[S,s]aprotroph",(ps@tax_table[,1]),value = TRUE) %>% 
  grep(pattern="[M,m]ycorrhizal",x=.,value = TRUE, invert = TRUE) %>% 
  unique()

# identify "pathogens"
pathogen_guilds <- 
  grep("[P,p]athogen|[P,p]arasite",(ps@tax_table[,1]),value = TRUE) %>% 
  grep(pattern="[M,m]ycorrhizal",x=.,value = TRUE, invert = TRUE) %>% 
  unique()

# subset taxa to only mutualists; get row sums; this will be proportion of mutualists
# in each sample

mutualist_proportions <- 
  ps %>% transform_sample_counts(function(x){x/sum(x)}) %>% 
  subset_taxa(Guild %in% mutualist_guilds) %>% 
  sample_sums()

# build data frame for modeling
mutualism_df <- 
  microbiome::meta(ps) %>% 
  mutate(proportion_mutualist = mutualist_proportions)

# subset taxa to only saprotrophs; get row sums; this will be proportion of mutualists
# in each sample

saprotroph_proportions <- 
  ps %>% transform_sample_counts(function(x){x/sum(x)}) %>% 
  subset_taxa(Guild %in% saprotroph_guilds) %>% 
  sample_sums()

# build data frame for modeling
saprotroph_df <- 
  microbiome::meta(ps) %>% 
  mutate(proportion_saprotroph = saprotroph_proportions)

# subset taxa to only pathogens; get row sums; this will be proportion of mutualists
# in each sample

pathogen_proportions <- 
  ps %>% transform_sample_counts(function(x){x/sum(x)}) %>% 
  subset_taxa(Guild %in% pathogen_guilds) %>% 
  sample_sums()

# build data frame for modeling
pathogen_df <- 
  microbiome::meta(ps) %>% 
  mutate(proportion_pathogen = pathogen_proportions)
pathogen_df

guild_df <- 
full_join(pathogen_df,mutualism_df) %>% 
  full_join(saprotroph_df)

## Fungal traits ####

# only got like 4 hits total...abandon this database for downstream work!

# # download traits metadata
# traits_meta <- read_csv("https://github.com/traitecoevo/fungaltraits/releases/download/v0.0.3/funtothefun.csv")
# 
# # download FungalTraits database
# traits_db <- fungaltraits::fungal_traits()
# names(traits_db$species)
# # match taxa at genus level
# genera <- ps@tax_table[,6] %>% str_remove("^g__")
# species <- ps@tax_table[,7] %>% str_remove("^s__")
# fungal_traits <- 
#   data.frame(Genus=genera) %>% 
#   mutate(species=paste(Genus,species,sep="_")) %>% 
#   left_join(traits_db,by=c("species","Genus"),multiple='all')
# View(fungal_traits)
# # need to condense/remove multiple matches
# fungal_traits %>% 
#   dplyr::filter(species != "NA_NA")
# 
# # remove traits not associated with biochem functional potential
# traits_to_ignore <- c(
#   "redChannel_mean","redChannel_sd","RNAHelicase_count","RNApolymerase_count","spore_length",
#   "spore_size","spore_width","sporocarp_chitin","sporocarp_N","sporocarp_protein","sporocarp_resp",           
#   "taxonomic_level_fg","tissue_c","tissue_cn","tissue_cp","tissue_n","tissue_np","tissue_p","total_genes",
#   "trehalase_count","latitude","map","greenChannel_mean","greenChannel_sd","heatShockProtein_count",
#   "extension_rate","fruiting_body_size","mat","longitude","melanin_content","melanin_count",
#   "coldShockProtein_count","dsDNA","blueChannel_mean","blueChannel_sd","ifungorum_number",
#   "sterol_type","studyName","substrate","trait_fg","trophic_mode_fg",'notes_fg',"source_funguild_fg",
#   "growth_form_fg","guild_fg","higher_clade","culture_media","culture_notes","elevation","em_expl",
#   "em_text","colour_mean","confidence_fg","ascoma_development","ascoma_type","ascus_dehiscence",
#   "uuid","obj_id","speciesMatched"
# )
# 
# # group by species; summarize to find mean values with na.omit=TRUE
# summarized_traits <- 
#   fungal_traits %>% 
#   dplyr::select(-all_of(traits_to_ignore)) %>% 
#   dplyr::group_by(species) %>% 
#   summarize(across(where(is.numeric),function(x){mean(x,na.rm=TRUE)}))
# 
# names(summarized_traits)
# 
# # join traits with tax_table species 
# 
# traits <- 
#   data.frame(Genus=genera) %>% 
#   mutate(species=paste(Genus,species,sep="_")) %>% 
#   left_join(summarized_traits,by=c("species"))
# View(traits)




# Models with watershed predictors ####


## Prepare responses ####

### Alpha-div ####
meta <- 
diversity %>%
  mutate(sampleid=row.names(.)) %>% 
  select(sampleid,Observed,Shannon) %>% # richness, Shannon div
  full_join(meta) 


# guilds/diversity
meta <- 
guild_df %>% 
  select(sampleid,starts_with("proportion_")) %>% 
  full_join(meta)


### Beta-div ####




### Gamma-div ####
# proportion of total taxa found locally (at each site)
n.taxa <- 
ps %>% 
  subset_taxa(taxa_sums(ps)>0) %>% 
  ntaxa()

meta$richness_proportion <- meta$Observed / n.taxa


### Rare/Unique taxa ####

# proportion of taxa below some relative abundance threshold
ra_threshold <- 0.001 # 0.1% total relative abundance threshold. Less than this is "rare"

n.taxa.total <- 
ps %>% 
  subset_taxa(taxa_sums(ps) > 0) %>% 
  taxa_sums() %>% sum()

taxa.proportions <- 
ps %>% 
  subset_taxa(taxa_sums(ps) > 0) %>% 
  taxa_sums() / n.taxa.total

(taxa.proportions < ra_threshold) %>% sum

rare_asvs <- 
  data.frame(
    asv = names(taxa.proportions),
    rare = (taxa.proportions < ra_threshold)
  ) %>% 
  dplyr::filter(rare) %>% 
  pluck("asv")

# proportion of total reads in each sample from "rare taxa"
sample.sums <- sample_sums(ps)

rare.asv.sums <- 
ps %>% 
  subset_taxa(taxa_names(ps) %in% rare_asvs) %>% 
  subset_taxa(taxa_sums(ps) > 0) %>% 
  sample_sums()

proportion_rare_taxa <- rare.asv.sums / sample.sums


meta <- 
data.frame(sampleid = names(proportion_rare_taxa),
           proportion_rare_taxa = proportion_rare_taxa) %>% 
  full_join(meta)

# meta %>% View

# proportion of taxa unique to a given site
pa <- 
ps %>% 
  otu_table() %>% 
  as("matrix") %>% 
  vegan::decostand(method = 'pa')

# find asvs that are unique to a single location
unique_asvs <- 
data.frame(unique = colSums(pa) < 2) %>% 
  mutate(asv = row.names(.)) %>%
  dplyr::filter(unique) %>% 
  pluck("asv")
# get total reads from those asvs for each sample
unique_sums <- 
ps %>% 
  subset_taxa(taxa_names(ps) %in% unique_asvs) %>% 
  sample_sums()
meta <- 
data.frame(
  proportion_unique_taxa = (unique_sums / sample.sums),
  sampleid = names((unique_sums / sample.sums))
) %>% 
  full_join(meta)


## Prepare predictors ####

### Watershed land use ####

read_csv("./watersheds/tables/nlcd_landcover.csv")


# read saved results
dflc = read.csv(file.path("./watersheds/tables", "nlcd_landcover.csv") ) %>% 
  mutate(location = site_ID)
dfid = read.csv(file.path("./watersheds/tables", "nlcd_impervious_desc.csv") ) %>% 
  mutate(location = site_ID)

# add land cover
meta <- 
dflc %>% 
  janitor::clean_names() %>% 
  dplyr::select(-site_id) %>% 
  full_join(meta)
meta <- 
dfid %>% 
  janitor::clean_names() %>% 
  dplyr::select(-site_id,-area_km2) %>% 
  full_join(meta)

# land use totals 
rsum.lc = dflc %>% select(-site_ID, -location, -area_km2, Developed..Open.Space, Developed..Low.Intensity,Developed..Medium.Intensity, Developed..High.Intensity,
) %>% rowSums()
total.developed <- 
  dflc %>% 
  select(-site_ID, -location, -area_km2, Developed..Open.Space, Developed..Low.Intensity,Developed..Medium.Intensity, Developed..High.Intensity) %>% 
  select(starts_with("Developed")) %>% 
  rowSums()
total.impervious <- dfid %>% select(-c(area_km2, site_ID, Unclassified,location)) %>%
  rowSums()

meta <- 
data.frame(location = seq_along(total.developed),
           total_developed_km2 = total.developed,
           total_impervious_km2 = total.impervious) %>% 
  full_join(meta)

# normalize land cover to proportions of total watershed area at each location
names(meta)

meta <- 
meta %>% 
  mutate(proportion_developed = total_developed_km2 / area_km2,
         proportion_impervious = total_impervious_km2 / area_km2,
         proportion_developed_high = developed_high_intensity / area_km2,
         proportion_developed_med = developed_medium_intensity / area_km2,
         proportion_developed_low = developed_low_intensity / area_km2,
         proportion_developed_openspace = developed_open_space / area_km2,
         proportion_cultivated_crops = cultivated_crops / area_km2)

# add land-use development categorical variable, selected by kmeans
set.seed(666)
km <- meta$proportion_developed %>% kmeans(2)
meta$developed_lgl <- (km$cluster == 2)

# find categorical groupings for development levels and impermeable surface area
set.seed(666)
km_developed_high <- meta$proportion_developed_high %>% kmeans(3)
km_impervious <- meta$proportion_impervious %>% kmeans(3)

meta <- 
data.frame(
  high_develop_cat = km_developed_high$cluster,
  impervious_cat = km_impervious$cluster
) %>% 
  mutate(high_develop_cat = case_when(high_develop_cat == 3 ~ "low",
                                      high_develop_cat == 1 ~ "med",
                                      high_develop_cat == 2 ~ "high"),
         impervious_cat = case_when(impervious_cat == 3 ~ "low",
                                    impervious_cat == 1 ~ "med",
                                    impervious_cat == 2 ~ "high"),
         across(everything(),function(x){factor(x,levels=c("low","med","high"))}),
         sampleid = meta$sampleid) %>% 
  full_join(meta)

### Water chemistry ####

# This was a no-go, with all contaminants testing at "too low to detect"



## Build models ####
meta %>% names

### gamma-div ####
gamma_div_mod <- 
glm(data = meta,
    formula = richness_proportion ~ proportion_developed + proportion_impervious + proportion_developed_high +
      proportion_developed_med + proportion_developed_low + proportion_developed_openspace + proportion_cultivated_crops)
gamma_div_step <- MASS::stepAIC(gamma_div_mod)
gamma_div_step$formula
gamma_div_mod <- 
  glm(data = meta,
      formula = gamma_div_step$formula)

summary(gamma_div_mod)
report::report(gamma_div_mod)

### shannon-div ####
shannon_div_mod <- 
glm(data = meta,
    formula = Shannon ~ proportion_developed + proportion_impervious + proportion_developed_high +
      proportion_developed_med + proportion_developed_low + proportion_developed_openspace + proportion_cultivated_crops)
shannon_div_step <- MASS::stepAIC(shannon_div_mod)
shannon_div_step$formula
shannon_div_mod <- 
  glm(data = meta,
      formula = shannon_div_step$formula)

summary(shannon_div_mod)

### rare-taxa ####
raretaxa_mod <- 
  glm(data = meta,
      formula = proportion_rare_taxa ~ proportion_developed + proportion_impervious + proportion_developed_high +
        proportion_developed_med + proportion_developed_low + proportion_developed_openspace + proportion_cultivated_crops)
raretaxa_step <- MASS::stepAIC(raretaxa_mod)
raretaxa_step$formula
raretaxa_mod <- 
  glm(data = meta,
      formula = raretaxa_step$formula)

summary(raretaxa_mod)

### unique-taxa ####
uniquetaxa_mod <- 
  glm(data = meta,
      formula = proportion_rare_taxa ~ proportion_developed + proportion_impervious + proportion_developed_high +
        proportion_developed_med + proportion_developed_low + proportion_developed_openspace + proportion_cultivated_crops)
uniquetaxa_step <- MASS::stepAIC(uniquetaxa_mod)
uniquetaxa_step$formula
uniquetaxa_mod <- 
  glm(data = meta,
      formula = uniquetaxa_step$formula)

summary(uniquetaxa_mod)


### NMDS / PermANOVA ####

# build new ps object with updated metadata
met <- sample_data(meta)
sample_names(met) <- meta$sampleid
sample_names(ps)

ps_meta <- phyloseq(otu_table(ps),
                    tax_table(ps),
                    sample_data(met),
                    phy_tree(ps))
ps_meta@sam_data %>% names
# ordinate with transformed data
ord <- ps_meta %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate(distance = 'unifrac')


scores <- vegan::rda(otu_table(ps_meta) %>% as("matrix"))
points1 <- scores$CA$u[,1]
points2 <- scores$CA$u[,2]
plot(ord$rproj[,1],ord$rproj[,2])
data.frame(DCA1=ord$rproj[,1],
           DCA2=ord$rproj[,2],
           proportion_developed=meta$proportion_developed,
           developed_category=meta$developed_lgl) %>% 
  ggplot(aes(x=DCA1,y=DCA2)) +
  geom_point(mapping=aes(color=proportion_developed,x=DCA1,y=DCA2),inherit.aes = FALSE) +
  # stat_ellipse(mapping=aes(color=developed_category,x=DCA1,y=DCA2),inherit.aes = FALSE) +
  labs(x="DCA1 [28.6%]",y="DCA2 [27.5%]",color="Impermeable surface\nproportion") +
  # scale_color_manual(values = pal.discrete[c(2,16)]) +
  theme_minimal() +
  theme(panel.grid = element_blank())

ggsave("./output/figs/DCA_Plot_part_1.png",dpi=300,height = 6,width = 6)
data.frame(DCA1=ord$rproj[,1],
           DCA2=ord$rproj[,2],
           proportion_developed=meta$proportion_developed,
           developed_category=meta$developed_lgl) %>% 
  ggplot(aes(x=DCA1,y=DCA2)) +
  # geom_point(mapping=aes(color=proportion_developed,x=DCA1,y=DCA2),inherit.aes = FALSE) +
  stat_ellipse(mapping=aes(color=developed_category,x=DCA1,y=DCA2),inherit.aes = FALSE) +
  labs(x="DCA1 [28.6%]",y="DCA2 [27.5%]",color="Highly developed") +
  scale_color_manual(values = pal.discrete[c(2,16)]) +
  theme_minimal() +
  theme(panel.grid = element_blank())
ggsave("./output/figs/DCA_Plot_part_2.png",dpi=300,height = 6,width = 6)

plot_ordination(ps_meta,ord,color = "developed_lgl") +
  stat_ellipse() +
  theme_minimal() +
  labs(color="Highly developed") 
  # scale_color_manual(values = pal.discrete[c(2,16)])
ggsave("./output/figs/DCA_Plot_land-use.png",dpi=300,width = 4,height = 4)

permanova_mod <- 
vegan::adonis2(data=meta,
               formula = otu_table(ps_meta) ~ meta$proportion_developed + meta$proportion_developed_low + meta$proportion_cultivated_crops)
as.data.frame(permanova_mod)







### guild proportions ####
mutualist_mod <- 
glm(data = meta,
    formula = proportion_mutualist ~ proportion_developed + proportion_impervious + proportion_developed_high +
      proportion_developed_med + proportion_developed_low + proportion_developed_openspace + proportion_cultivated_crops)
mutualist_step <- MASS::stepAIC(mutualist_mod)
mutualist_step$formula
mutualist_mod <- 
  glm(data = meta,
      formula = mutualist_step$formula)
summary(mutualist_mod)

pathogen_mod <- 
  glm(data = meta,
      formula = proportion_pathogen ~ proportion_developed + proportion_impervious + proportion_developed_high +
        proportion_developed_med + proportion_developed_low + proportion_developed_openspace + proportion_cultivated_crops)
pathogen_step <- MASS::stepAIC(pathogen_mod)
pathogen_step$formula
pathogen_mod <- 
  glm(data = meta,
      formula = pathogen_step$formula)
summary(pathogen_mod)

### Save model objects ####
mod_objects <- c("gamma_div_mod",
                 "shannon_div_mod",
                 "uniquetaxa_mod",
                 "raretaxa_mod",
                 "permanova_mod",
                 "mutualist_mod",
                 "pathogen_mod")

for(i in mod_objects){
  fn <- paste0("./output/models/",i,".RDS")
  saveRDS(object = get(i),file = fn)
}


# Differential abundance ####
# quick and dirty diffabund analysis
# This should use categorical Urban/Non-Urban predictor!
    # Define that based on watershed data #

ps_meta@sam_data %>% names

da_analysis_hi_devel <- differentialTest(formula = ~ high_develop_cat, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", boot = FALSE,
                                         data = ps_meta,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)

da_analysis_hi_devel$significant_taxa

bbdml_hi_devel <- multi_bbdml(da_analysis_hi_devel,ps_meta,"high_develop_cat","high_develop_cat")

plot_multi_bbdml(bbdml_hi_devel,color = "high_develop_cat")
saveRDS(bbdml_plot_1,"./output/figs/diffabund_plot_high_development_categorical.RDS")

# check this taxon for more current names
taxon <- 
otu_to_taxonomy(da_analysis_hi_devel$significant_taxa,ps) %>% 
  str_split("g__") %>% map_chr(2) %>% str_replace("_s__"," ")
mb_db <- mycobank::get_mycobank_db()
mycobank::get_mycobank_synonyms(taxon,mb_db)$current_name


da_analysis_imperm <- differentialTest(formula = ~ impervious_cat, #abundance
                                         phi.formula = ~ 1, #dispersion
                                         formula_null = ~ 1, #mean
                                         phi.formula_null = ~ 1,
                                         test = "Wald", boot = FALSE,
                                         data = ps_meta,
                                         fdr_cutoff = 0.05,
                                         full_output = TRUE)

da_analysis_imperm$significant_taxa

bbdml_imperm <- multi_bbdml(da_analysis_hi_devel,ps_meta,"impervious_cat","impervious_cat")
plot_multi_bbdml(bbdml_imperm,color = "impervious_cat")
saveRDS(bbdml_plot_1,"./output/figs/diffabund_plot_impervious_surface_categorical.RDS")
# https://www.coelomycetes.org/pleosporales/neopyrenochaetaceae/neopyrenochaeta/neopyrenochaeta-annellidica-.html

# glm along gradient for significant taxa

da_taxa <- otu_to_taxonomy(da_analysis_imperm$significant_taxa,ps_meta) %>% str_split("_") %>% unlist()
da_taxa_name <- paste0(da_taxa[15]," ",da_taxa[18])
meta[,da_taxa_name] <- 
  ps_meta %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  subset_taxa(taxa_names(ps_meta) %in% da_analysis_imperm$significant_taxa) %>% 
  sample_sums()

N.annelidica_mod <- 
glm(data = meta,
    formula = `Neopyrenochaeta annellidica` ~ proportion_developed + proportion_impervious + proportion_developed_high +
      proportion_developed_med + proportion_developed_low + proportion_developed_openspace + proportion_cultivated_crops)
N.annelidica_step <- MASS::stepAIC(N.annelidica_mod)
N.annelidica_mod <- glm(data=meta,
                        formula = N.annelidica_step$formula)
summary(N.annelidica_mod)

saveRDS(N.annelidica_mod,"./output/models/N.annelidica_model.RDS")


# save metadata object for further plotting
saveRDS(meta,"./output/cleaned_full_metadata.RDS")

# count all downstream taxa from each successive point
# taxa from top site not found in any site below
# taxa from 1 & 2 not found in 3-10, etc

site.unique.taxa <- list()
for(i in 10:1){
  x <- ps %>% 
    subset_samples(location == i)
  x <- x %>% 
    subset_taxa(taxa_sums(x) > 0)
  site.unique.taxa[[paste0("location_",as.character(i))]] <- taxa_names(x)
}
site.unique.taxa$location_1

grouplist <- list()
for( i in 10:1){
  groupname <- paste0(i," - 1")
  grouplist[[groupname]]  <- site.unique.taxa[1:i] %>% unlist %>% unique
  # purrr::reduce(unique)
}
data.frame(downstream_position = 1:10,
           total_taxa_downstream = 
             grouplist %>% map_dbl(length)) %>% 
  mutate(total_taxa_upstream = max(total_taxa_downstream) - total_taxa_downstream) %>% 
  pivot_longer(starts_with("total_"),names_to = "position", values_to = "n_unique_taxa",
               names_prefix = "total_taxa_") %>% 
  mutate(downstream_position = factor(downstream_position, levels = as.character(1:10))) %>% 
  ggplot(aes(x=downstream_position,y=n_unique_taxa/max(n_unique_taxa),fill=position)) +
  geom_col() +
  theme_bw() +
  labs(x="Location",y="Proportion of taxa unique\n to cumulative watershed") +
  scale_fill_viridis_d(end=.8)
ggsave("./output/figs/n_unique_taxa_upstream_and_downstream.png")
