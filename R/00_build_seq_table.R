# -----------------------------------------------------------------------------#
# Utah county urban stream endophytes ITS DADA2 pipeline
# Processing raw amplicon reads (these have already been passed through itsxpress)
# Author: Geoffrey Zahn
# Software versions:  R v 4.2.2
#                     tidyverse v 1.3.2
#                     dada2 v 1.24.0
#                     decontam v 1.16.0
#                     phyloseq v 1.40.0
#                     Biostrings v 2.64.0
#                     patchwork v 1.1.1
#                     readxl v 1.4.1
#                     janitor::clean_names() v 2.1.0
# -----------------------------------------------------------------------------#

#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal  # 
#################################################################################

# PACKAGES, SCRIPTS, AND SETUP ####

# why each package (put in onboarding document)
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")

set.seed(666) # "random" seed for reproducibility


# PARSE FILE PATHS ####

# File parsing - 

path <- "./data" # CHANGE to the directory containing your demultiplexed fastq files when using your own data
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present


# Here, the pattern is set to match the forward reads tagged as "...pass_1.fastq.gz" in their filenames; 
# "...pass_2.fastq.gz" for reverse reads

# Your data may differ, using "F" and "R" in the filenames, or something similar..
# Be sure to change that pattern to fit your own files when using your own data
fns <- sort(list.files(file.path(path,"raw"), full.names = TRUE, pattern = "ITS1.fastq.gz")) # make pattern match your FWD reads

sample.names <- basename(fns) %>% str_remove(".raw_1.ITS1.fastq.gz")

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
# you can select any number of files here...
# as-is, this just shows the Fwd and Rev quality profiles for the 1st and 2nd files
p1 <- plotQualityProfile(fns[1]) + ggtitle("Example forward reads")

# display and save the plots
p1
ggsave("./output/figs/unfiltered_quality_plots.png",dpi=500,height = 6,width = 6)

# FILTER AND TRIM ####

# here, we decide what filenames we want our filtered reads to be called
# in this case, we use the base name of the sample and save the filtered files into the "filtered" subdirectory
filts_f <- file.path(path, "filtered", paste0(sample.names, "_FWD_filt.fastq.gz"))

# this is the actual qualit control step
# These values are informed by our quality plots
out <- filterAndTrim(fns, filts_f, # input and output file names as denoted above
                     maxN=0, # uncalled bases are currently not supported in dada2
                     maxEE=2, # refers to the maximum expected errors allowed
                     truncQ=2, # special value denoting "end of good quality sequence" (optional)
                     rm.phix=TRUE, # automatically remove PhiX spike-in reads from sequencing center
                     compress=TRUE, # compress output files with gzip
                     multithread=4) # On Windows set multithread=FALSE

# save the filtration info in case we need it later
saveRDS(out, "./output/trackreads.RDS")


# In case some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
sample.names <- basename(filts_f) %>% str_remove("_FWD_filt.fastq.gz")


# sanity check  comparison of before and after filtration
p3 <- plotQualityProfile(fns[1]) + ggtitle("Unfiltered")
p4 <- plotQualityProfile(filts_f[1])+ ggtitle("Filtered")
p3 / p4
ggsave("./output/figs/filtered_quality_comparison.png",dpi=300,height = 6,width = 6)

# LEARN ERROR RATES ####

# learn errors
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20,verbose = 2) # set multithread = FALSE on Windows
saveRDS(errF,"./output/errF.RDS")

# sanity check for error model
# explain what to look for in the plot! 
plotErrors(errF, nominalQ=FALSE)
ggsave("./output/figs/error_model.png",dpi=200,height = 6,width = 6)

# dereplication
derepF <- derepFastq(filts_f, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepF) <- 
  names(derepF) %>% 
  str_remove("_FWD_filt.fastq.gz")

# SAMPLE INFERRENCE ####
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo") # set multithread = FALSE on Windows
saveRDS(dadaFs,"output/dadaFs.RDS")


# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(dadaFs)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"./output/seqtab.nochim.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF","nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "./output/ITS_read_counts_at_each_step.csv", row.names = TRUE)

############# PICK UP HERE ~~~~~~~~~~~~~~~~~~~~~~~~~

# make metadata, for now
meta <- 
data.frame(sampleid=row.names(track),
           control=row.names(track) %>% grepl(pattern="CTL"))





# IMPORT METADATA ####

# import and clean
# meta <- read_xlsx("./data/SI_Indo_Metadata_CorrectGPS.xlsx") %>%
#   janitor::clean_names() %>%
#   dplyr::rename("east_west" = "east_or_west_of_wallace_line") %>%
#   mutate(gps=case_when(gps == "N/A" ~ "NA NA",
#                        TRUE ~ gps),
#          lat=gps %>% str_split(" ") %>% map_chr(1) %>% as.numeric(),
#          lon=gps %>% str_split(" ") %>% map_chr(2) %>% as.numeric(),
#          blank = grepl("Blank",ignore.case = TRUE, sample)) %>%
#   select(sample, location, east_west, lat, lon, blank)
# 
# # subset to match seq table sample names
# meta <- meta[meta$sample %in% sample.names, ]

# reorder metadata to match seqtab
# df <- data.frame(seqtab_rows=row.names(seqtab.nochim),
#                  sample=row.names(seqtab.nochim))
# df2 <- left_join(meta,df,by="sample")
# row.names(df2) <- df2$sample
# meta <- df2[row.names(seqtab.nochim),]
# row.names(meta) <- meta$sample
# identical(row.names(meta),row.names(seqtab.nochim))
# 

# Remove all seqs with fewer than 100 nucleotides (if any) ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./output/seqtab.nochim.clean.RDS")

# Find and remove contaminants ####
contams = isContaminant(seqtab.nochim, neg = meta$control, normalize = TRUE)

table(contams$contaminant) # how many taxa are contaminants?
write.csv(contams, file = "./output/likely_contaminants.csv", row.names = TRUE)

# remove contaminant sequences and control samples from both tables, respectively ####
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$control == FALSE,]
meta = meta[meta$control == FALSE,]
dim(seqtab.nochim)


# reload point
saveRDS(seqtab.nochim,"./output/seqtab_for_taxonomy.RDS")

