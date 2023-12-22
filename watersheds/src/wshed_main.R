#
# Workflow for watershed delineation
#
source("./wshed_functions.R")

#=========================================================
# Main - Watershed delineation / NLCD stats
#=========================================================
#
# Required data: 
#   - DEM (10-meter 3DEP used), 
#   - CSV file with GPS locations
#   - NLCD landcover and impervious descriptor layers (2021 used) to extent of DEM
# optional: 
#   - HUC8 watershed for general area
#
## Run paths
# set up directory to save/work
outpath = "."
dempath = "."
nlcdpath <- "."
csvpath <- ".s/sample-sites.csv" # original csv file
shpath = "./sample_sites.shp" # created

## Run functions
# create vect object for points
geospat_csv(csvpath, write.path=".")
# several watershed related functions for delineation with DEM and NLCD layer masking
wbt_watersheds(dempath, shpath, outpath)
# create data table with final results
extract_lc_table(outpath)

# note: several new files and folders will be saved in outpath location

#=========================================================
# Plotting 
#=========================================================
#
# # # # # # # # # # # # # # # # # # # # # # # 
#
# a few simple plots with these results
#
library(gridExtra);library(ggplot2);library(reshape2)
# read saved results
dflc = read.csv(file.path(outpath, "tables", "nlcd_landcover.csv") )
dfid = read.csv(file.path(outpath, "tables", "nlcd_impervious_desc.csv") )

# row sums 
rsum.lc = dflc %>% select(site_ID, Developed..Open.Space, Developed..Low.Intensity,Developed..Medium.Intensity, Developed..High.Intensity,
) %>% rowSums()
rsum.id = dfid %>% select(-c(area_km2, site_ID, Unclassified)) %>%
  rowSums()

# Figure comparison
p1= dflc %>% select(site_ID, Developed..Open.Space, Developed..Low.Intensity,Developed..Medium.Intensity, Developed..High.Intensity) %>% 
  mutate(Total.Developed = rsum.lc) %>%  melt(.,id.vars="site_ID") %>% mutate(Classification = variable) %>%
  ggplot(aes(x=site_ID,y=value,colour=Classification, linetype=Classification)) + geom_point() +
  geom_line() + theme_classic() + ylab(expression("Developed Land Cover Area ( km"^2~")") ) +
  scale_x_continuous(breaks=seq(1,10)) + xlab("") + 
  scale_linetype_manual(values=c(rep("dashed",4),"solid"))

p2 = dfid %>% select(-c(area_km2,Unclassified)) %>% mutate(Total.impervious = rsum.id) %>% 
  melt(.,id.vars="site_ID") %>% mutate(Classification = variable) %>% 
  ggplot(aes(x=site_ID,y=value,colour=Classification, linetype=Classification)) + geom_point() +
  geom_line() + theme_classic() + scale_color_hue(l=50,c=80) + 
  ylab(expression("Impervious Surface Area ( km"^2~")") ) +
  scale_x_continuous(breaks=seq(1,10)) + xlab("Site ID") + 
  scale_linetype_manual(values=c(rep("dashed",5),"solid"))

# plot impervious
p2 + ggtitle("Impervious Surface Area by Site Watershed (2021)")

# plot together
grid.arrange(p1,p2,nrow=2)

# ## Save
# dir.create(file.path(outpath,"figures"))
# # save 1
# png(file.path(outpath,"figures","figure1.png"), width=9,height=6, units='in', res=300)
# p2 + ggtitle("Impervious Surface Area by Site Watershed (2021)")
# dev.off()
# # save 2
# png(file.path(outpath,"figures","figure2.png"), width=8.5,height=5.5, units='in', res=300)
# grid.arrange(p1,p2,nrow=2)
# dev.off()

# map
