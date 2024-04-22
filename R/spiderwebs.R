library(readxl)
library(tidyverse)
library(vegan)
library(ggmap)
source("./googlemap_styling.R")

mapstyle <- rjson::fromJSON(file = "./mapstyle3.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

set.seed(666)

# Load google maps API key from .Renviron 
ggmap::register_google(key = Sys.getenv("APIKEY")) # Key kept private


dat <- read_xlsx("./fungi data with species.xlsx",sheet = 2) %>% 
  filter(!is.na(LAT)) %>% 
  mutate(set = case_when(SITE %in% c("A","B") ~ "AB",
                         SITE %in% c("C","D") ~ "CD",
                         SITE %in% c("E","F") ~ "EF"),
         status = case_when(SITE %in% c("A","C","E") ~ "Burned",
                            TRUE ~ "Unburned"))

dat[is.na(dat)] <- 0

area <- 
  ggmap::get_googlemap(center = c(lon = mean(dat$LON,na.rm = TRUE), 
                                  lat = mean(dat$LAT,na.rm=TRUE)),
                       zoom = 11,
                       scale = 2,
                       style=mapstyle)
ggmap(area) +
  geom_point(data=dat,aes(x=LON,y=LAT),color='darkorange',size=.5,alpha=.75) 
ggsave("./location_map.png",dpi=300,height = 8,width = 8)


# make matrix of species presence
mat <- 
dat %>% 
  select(-c(ID,SITE,NUMBER,LAT,LON,DATE,set,status))


mat <- mat[which(rowSums(mat) != 0),]
dat <- dat[which(rowSums(mat) != 0),]



nmds <- metaMDS(as.matrix(mat))
mds1 <- nmds$points[,1]
mds2 <- nmds$points[,2]

data.frame(location = dat$set,status = dat$status,
           MDS1 = mds1, MDS2 = mds2) %>% 
  ggplot(aes(x=mds1,y=mds2,color=status)) +
  geom_point(size=3) +
  stat_ellipse() +
  scale_color_manual(values = c("#fc2c03","#38a336")) +
  theme_minimal() +
  theme(legend.text = element_text(face='bold',size=20),
        legend.title  = element_text(face='bold',size=24))
ggsave("./nmdsplot.png",dpi=300,height = 8,width = 8)
adonis2(mat ~ dat$status)
