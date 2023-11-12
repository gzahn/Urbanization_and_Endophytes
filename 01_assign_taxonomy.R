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
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, combine sequence table and metadata                     #
#################################################################################


library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(patchwork); packageVersion("patchwork")


# Load seq table
seqtab <- readRDS("./output/seqtab_for_taxonomy.RDS")

colnames(seqtab)[1]
seqtab[1,]
taxa <- assignTaxonomy(colnames(seqtab)[1], "./taxonomy/sh_general_release_dynamic_s_all_16.10.2022.fasta.gz", multithread=4)


# assign taxonomy
# Use Unite_All
taxa <- assignTaxonomy(seqtab, "./taxonomy/sh_general_release_dynamic_s_all_16.10.2022.fasta.gz", multithread=4)

# Save taxonomy file
saveRDS(taxa, file = "./output/RDP_Taxonomy_from_dada2.RDS")


# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(meta)
sample_names(met) <- meta$sample

ps <- phyloseq(otu,met,tax)
saveRDS(ps,"./output/ps_not-cleaned.RDS")
ps
