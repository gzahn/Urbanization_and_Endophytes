# Build phyloseq object

library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(FUNGuildR)

# read otu table
seqtab.nochim <- readRDS("./output/seqtab_for_taxonomy.RDS")

# read taxonomy
taxa <- readRDS("./output/RDP_Taxonomy_from_dada2.RDS")



# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
# met <- sample_data(meta)

meta <- 
data.frame(sampleid=row.names(seqtab.nochim),
           location=row.names(seqtab.nochim) %>% 
             str_remove("JM_") %>% 
             str_remove("A") %>% 
             str_remove("B") %>% 
             str_remove("C"))

met <- sample_data(meta)
sample_names(met) <- meta$sampleid

ps <- phyloseq(otu,met,tax)
saveRDS(ps,"./output/ps_not-cleaned_no-metadata.RDS")
ps


