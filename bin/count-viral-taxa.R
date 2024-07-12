#!/usr/bin/env Rscript

library(optparse)
library(tidyverse)

# Set arguments
option_list = list(
  make_option(c("--reads"), type="character", default=NULL,
              help="Path to TSV file of human-infecting virus read data."),
  make_option(c("--taxa"), type="character", default=NULL,
              help="Path to database of viral taxonomic relationships"),
  make_option(c("--output"), type="character", default=NULL,
              help="Path to output file.")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Set input and output paths
read_db_path <- opt$reads
taxa_db_path <- opt$taxa
out_path   <- opt$output

#=====================#
# AUXILIARY FUNCTIONS #
#=====================#

count_children <- function(taxid_tree){
    # Count number of direct children of each node in a taxid tree
    n_children <- taxid_tree %>% group_by(parent_taxid) %>%
        count(name="n_children") %>% rename(taxid = parent_taxid)
    db_out <- taxid_tree %>% left_join(n_children, by="taxid") %>%
        mutate(n_children = replace_na(n_children, 0))
    return(db_out)
}

count_hits_direct <- function(read_db, taxid_tree, group_var){
    # Count number of reads directly assigned to each taxid in a tree
    taxids <- count_children(taxid_tree)
    direct_hits_setup <- read_db %>% group_by(taxid_best, .data[[group_var]]) %>%
        count(name="n_reads_direct", .drop=FALSE) %>%
        pivot_wider(names_from=any_of(group_var), values_from=n_reads_direct,
                    values_fill = 0)
    direct_hits <- taxids %>% left_join(direct_hits_setup, by=c("taxid" = "taxid_best")) %>%
        pivot_longer(cols=-(taxid:n_children), names_to=group_var,
                     values_to="n_reads_direct") %>%
        mutate(n_reads_direct = replace_na(n_reads_direct, 0))
    # Double pivot ensures all viruses have a count for all samples
    return(direct_hits %>% select(-n_children))
}

count_hits_indirect <- function(direct_hits, taxid_tree, group_var){
    # Count number of reads assigned to each taxid in a tree or its descendents
    iter <- 0
    taxids <- count_children(taxid_tree)
    n_children <- taxids %>% select(taxid, n_children)
    hits_unresolved <- direct_hits %>% left_join(n_children, by="taxid") %>%
        mutate(n_reads_clade = n_reads_direct, n_reads_children = 0)
    hits_resolved <- tibble()
    cat("\tIteration ", iter, ": ", nrow(hits_unresolved), " unresolved, ", nrow(hits_resolved), " resolved.\n", sep="")
    while (nrow(hits_unresolved) > 0){
        iter <- iter + 1
        # Identify leaf nodes
        hits_leaf <- hits_unresolved %>% filter(n_children == 0)
        # Calculate contributions to parent nodes
        hits_children <- hits_leaf %>% group_by(parent_taxid, .data[[group_var]]) %>%
            summarize(n_reads_children = sum(n_reads_clade), .groups = "drop") %>%
            rename(taxid = parent_taxid)
        # Add child reads to direct hits to calculate clade reads
        hits_partial <- hits_unresolved %>% filter(n_children > 0) %>%
            select(-n_children, -n_reads_children) %>%
            left_join(hits_children, by=c("taxid", group_var)) %>%
            mutate(n_reads_children = replace_na(n_reads_children, 0),
                   n_reads_clade = n_reads_clade + n_reads_children)
        # Add leaf nodes to resolved hits
        hits_resolved <- hits_resolved %>% bind_rows(hits_leaf) %>%
            select(-n_children, -n_reads_children)
        # Recalculate child counts after removing leaf nodes
        hits_n_children <- hits_partial %>% group_by(parent_taxid, .data[[group_var]]) %>%
            count(name="n_children") %>% rename(taxid=parent_taxid)
        hits_unresolved <- hits_partial %>%
            left_join(hits_n_children, by=c("taxid", group_var)) %>%
            mutate(n_children = replace_na(n_children, 0))
        cat("\tIteration ", iter, ": ", nrow(hits_unresolved), " unresolved, ", nrow(hits_resolved), " resolved.\n", sep="")
    }
    return(hits_resolved %>% arrange(taxid))
}

count_hits <- function(read_db, taxid_tree, group_var){
    # Count direct and indirect hits for each taxon in a set of human-viral reads
    cat("Counting direct hits...\n")
    hits_direct <- count_hits_direct(read_db, taxid_tree, group_var)
    cat("Done.\n")
    cat("Counting indirect hits...\n")
    hits_indirect <- count_hits_indirect(hits_direct, taxid_tree, group_var)
    cat("Done.\n")
    return(hits_indirect)
}

#============#
# RUN SCRIPT #
#============#

# Import data
read_db <- read_tsv(read_db_path, show_col_types = FALSE) %>% mutate(taxid_best = as.integer(taxid_best))
taxa_db <- read_tsv(taxa_db_path, show_col_types = FALSE) %>% mutate(taxid = as.integer(taxid))

if (nrow(read_db) > 0){
    # Count hits
    hits_db <- count_hits(read_db, taxa_db, "sample")
    # Filter out null rows
    hits_db_filtered <- filter(hits_db, n_reads_clade > 0)
} else {
    hits_db_filtered <- tibble(
        taxid = numeric(), name = character(), rank = character(), parent_taxid = numeric(),
        sample = character(), n_reads_direct = numeric(), n_reads_clade = numeric()
    )
}

# Write output
write_tsv(hits_db_filtered, out_path)
