## name: 1_btp_construct_ps
## purpose: to construct the phyloseq object combining the sequence 
## counts table, taxonomic annotations, and sample metadata
## author: Dr Aditya Bandla
## email: adityabandla@u.nus.edu

## libraries
suppressPackageStartupMessages({
  library(dada2)
  library(Hmisc)
  library(janitor)
  library(lubridate)
  library(magrittr)
  library(phyloseq)
  library(tidyverse)
})

## paths to directories
repo <- file.path("/Users/abandla/Desktop/2_research/2_manuscripts/2_2020_brunei_tropical_peat/")
data <- file.path(repo, "1_data")
constructs <- file.path(data, "3_constructs")
figures <- file.path(repo, "3_figures")

## export sequences
readRDS(file.path(data, "1_seqtab", "2020_btp_seqtab_nc.rds")) %>%
  getUniques(.) %>%
  uniquesToFasta(
    .,
    file.path(constructs, "2020_btp_asv_seqs.fa"),
    ids = paste0("ASV_", seq(length(getUniques(.))))
  )

## import sequence table
seqtab <- readRDS(file.path(data, "1_seqtab", "2020_btp_seqtab_nc.rds")) %>%
  set_colnames(paste0("ASV_", seq(length(getSequences(.))))) %>%
  set_rownames(gsub("\\-", "_", rownames(.)))

## import sample data
metadata <- read.csv(
  file.path(constructs, "1_metadata", "2020_btp_metadata.csv"),
  header = TRUE, row.names = 1, stringsAsFactors = TRUE
)

## import taxonomic data
taxa <- readRDS(file.path(data, "1_seqtab", "2020_btp_taxa_silva_v138.rds")) %>%
  set_rownames(paste0("ASV_", seq(length(getSequences(.)))))

## construct phyloseq object
btp_ps <- phyloseq(
  otu_table(seqtab, taxa_are_rows = FALSE),
  sample_data(metadata),
  tax_table(taxa)
) %>%
  subset_taxa(., !is.na(Kingdom)) %>%
  subset_taxa(., Order %nin% c("Chloroplast") & Family %nin% c("Mitochondria"))

## save phyloseq object
saveRDS(btp_ps, file.path(constructs, "2_phyloseq", "2020_btp_ps.rds"))

## subset phyloseq object to peat samples from the burnt 
## and intact psf study
btp_fire_ps <- btp_ps %>%
  subset_samples(psf_type != "NA") %>%
  filter_taxa(., function(x) sum(x > 0) >= 2, TRUE)

## save phyloseq object
saveRDS(btp_fire_ps, file.path(constructs, "2_phyloseq", "2020_btp_fire_ps.rds"))

## subset phyloseq object to peat samples from the 
## phytobiome study
btp_root_ps <- btp_ps %>%
  subset_samples(compartment != "NA") %>%
  filter_taxa(., function(x) sum(x > 0) >= 2, TRUE)

## inspect phyloseq object
btp_root_ps

## save phyloseq object
saveRDS(btp_root_ps, file.path(constructs, "2_phyloseq", "2020_btp_root_ps.rds"))

## libraries
suppressPackageStartupMessages({
  library(DECIPHER)
  library(phyloseq)
  library(tidyverse)
})

## import phyloseq object
btp_fire_ps <- readRDS("2020_btp_fire_ps.rds")

## import sequences and filter
readDNAStringSet("2020_btp_asv_seqs.fa") %>%
  .[names(.) %in% taxa_names(btp_fire_ps)] %>%
  writeXStringSet(., "2020_btp_fire_asv_seqs.fa")

## align sequences
readDNAStringSet("2020_btp_asv_seqs.fa") %>%
  .[names(.) %in% taxa_names(btp_fire_ps)] %>%
  AlignSeqs(., anchor = NA, processors = 24) %>%
  writeXStringSet(., "2020_btp_fire_asv_seqs_aligned.fa")

## build phylogenetic tree
## FastTree -nt 2020_btp_fire_asv_seqs_aligned.fa > 2020_btp_fire_asv_seqs_aligned.tre
