#!/usr/bin/env Rscript

library(jsonlite)
library(optparse)
library(tidyverse)

# Set arguments
option_list = list(
  make_option(c("--nodes"), type="character", default=NULL,
              help="Path to nodes DMP file giving taxonomy structure."),
  make_option(c("--names"), type="character", default=NULL,
              help="Path to names DMP file giving taxid-name mappings."),
  make_option(c("--hv"), type="character", default=NULL,
              help="Path to TSV containing human-infecting virus taxids."),
  make_option(c("--json"), type="character", default=NULL,
              help="Path to JSON file mapping viral taxids to descendent taxids from gimme-taxa.py."),
  make_option(c("--output"), type="character", default=NULL,
              help="Path to output file.")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Set input and output paths
nodes_path <- opt$nodes
names_path <- opt$names
hv_path    <- opt$hv
json_path  <- opt$json
out_path   <- opt$output

#=====================#
# AUXILIARY FUNCTIONS #
#=====================#

expand_taxids <- function(taxids_in, nodes){
  # Get all taxids descended from a list of input taxids
  taxids_out <- taxids_in
  taxids_new <- filter(nodes, parent_taxid %in% taxids_out, !taxid %in% taxids_out) %>%
      pull(taxid) %>% sort
  while (length(taxids_new) > 0){
      taxids_out <- c(taxids_out, taxids_new) %>% sort
      taxids_new <- filter(nodes, parent_taxid %in% taxids_out, 
                           !taxid %in% taxids_out) %>%
          pull(taxid) %>% sort
  }
  return(taxids_out)
}

#============#
# RUN SCRIPT #
#============#

# Import data
nodes <- read_tsv(nodes_path, col_names = FALSE) %>%
    transmute(taxid = X1, parent_taxid = X3, rank = X5)
names <- read_tsv(names_path, col_names = FALSE) %>%
    transmute(taxid = X1, name = X3, alt_name = X5, type = X7)
hv    <- read_tsv(hv_path, col_names <- c("taxid", "name"))
hv_extra <- read_json(json_path)

# Get viral taxid list
viral_taxids <- expand_taxids(c(10239, hv$taxid), nodes)

# Generate viral taxid DB
viral_db <- filter(nodes, taxid %in% viral_taxids) %>%
    left_join(names, by = "taxid") %>% select(-alt_name)
viral_db_scinames <- filter(viral_db, type == "scientific name") %>%
    group_by(taxid) %>% filter(row_number() == 1) %>% select(-type)
viral_db_other <- filter(viral_db, ! taxid %in% viral_db_scinames$taxid) %>%
    group_by(taxid) %>% filter(row_number() == 1) %>% select(-type)
viral_db_out <- bind_rows(viral_db_scinames, viral_db_other) %>%
    select(taxid, name, rank, parent_taxid) %>%
    arrange(taxid)

# Extract & process additional taxid from NCBI
taxids_extra <- lapply(names(hv_extra), function(taxid) 
                       tibble(ancestor_taxid=taxid, descendent_taxid = unlist(hv_extra[[taxid]]))) %>%
    bind_rows()
taxids_extra_filtered <- filter(taxids_extra, ! descendent_taxid %in% viral_db_out$taxid) %>%
    rename(taxid = ancestor_taxid) %>% mutate(taxid = as.integer(taxid)) %>%
    left_join(viral_db_out, by = "taxid") %>%
    group_by(descendent_taxid) %>% 
    filter(! taxid %in% parent_taxid) # Remove redundant ancestors, keeping only lowest-ranked
taxids_extra_single <- taxids_extra_filtered %>% group_by(descendent_taxid) %>% 
    filter(n() == 1) %>%
    select(taxid = descendent_taxid, parent_taxid = taxid) %>%
    mutate(name = NA, rank = NA)
# If still multiple possible ancestors, pick the lowest rank
ranks <- c("subspecies", "species", "subgenus", "genus", 
           "subfamily", "family", "suborder", "order", "class",
           "subphylum", "phylum", "kingdom", "superkingdom")
taxids_extra_multi <- taxids_extra_filtered %>% 
    filter(!descendent_taxid %in% taxids_extra_single$taxid) %>%
    mutate(rank_index = match(rank, ranks),
           rank_index = replace_na(rank_index, 0)) %>%
    group_by(descendent_taxid) %>%
    filter(rank_index == min(rank_index)) %>%
    filter(taxid == max(taxid)) %>%
    filter(n() == 1) %>%
    select(taxid = descendent_taxid, parent_taxid = taxid) %>%
    mutate(name = NA, rank = NA)
taxids_extra_out <- bind_rows(taxids_extra_single, taxids_extra_multi)

# Bind rows and write
taxids_all_out <- bind_rows(viral_db_out, taxids_extra_out) %>%
    arrange(taxid) %>%
    mutate(rank = replace_na(rank, "no rank"),
           name = replace_na(name, "unknown virus"))
write_tsv(taxids_all_out, out_path)
