// Further expand human-virus DB with additional taxa from NCBI
process FINALIZE_HUMAN_VIRUS_DB {
    input:
        path(human_virus_db)
        path(human_virus_descendents_json)
        path(taxonomy_nodes)
        path(taxonomy_names)
    output:
        path("human-virus-db.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript

        # Setup
        library(jsonlite)
        library(optparse)
        library(tidyverse)
        nodes_path <- "!{taxonomy_nodes}"
        names_path <- "!{taxonomy_names}"
        hv_path    <- "!{human_virus_db}"
        json_path  <- "!{human_virus_descendents_json}"
        out_path   <- "human-virus-db.tsv.gz"

        # Import data
        nodes <- read_tsv(nodes_path, col_names = FALSE) %>%
            transmute(taxid = X1, parent_taxid = X3, rank = X5)
        names <- read_tsv(names_path, col_names = FALSE) %>%
            transmute(taxid = X1, name = X3, alt_name = X5, type = X7)
        hv    <- read_tsv(hv_path, col_names <- c("taxid", "name"))
        hv_extra <- read_json(json_path)

        # Generate taxid DB from nodes and names
        hv_db <- filter(nodes, taxid %in% hv$taxid) %>%
            left_join(names, by="taxid") %>% select(-alt_name)
        hv_db_scinames <- filter(hv_db, type == "scientific name") %>%
            group_by(taxid) %>% filter(row_number() == 1) %>% select(-type)
        hv_db_other <- filter(hv_db, ! taxid %in% hv_db_scinames$taxid) %>%
            group_by(taxid) %>% filter(row_number() == 1) %>% select(-type)
        hv_db_out <- bind_rows(hv_db_scinames, hv_db_other) %>%
            select(taxid, name, rank, parent_taxid) %>% arrange(taxid)

        # Extract and process additional descendent taxids from NCBI
        ranks <- c("subspecies", "species", "subgenus", "genus", 
                   "subfamily", "family", "suborder", "order", "class",
                   "subphylum", "phylum", "kingdom", "superkingdom")
        taxids_extra <- lapply(names(hv_extra), function(taxid)
                                tibble(ancestor_taxid=taxid, descendent_taxit=unlist(hv_extra[[taxid]]))) %>%
            bind_rows()
        taxids_extra_filtered <- filter(taxids_extra, !descendent_taxid %in% hv_db_out$taxid) %>%
            rename(taxid = ancestor_taxid) %>% mutate(taxid = as.integer(taxid)) %>%
            left_join(hv_db_out, by = "taxid") %>%
            group_by(descendent_taxid) %>% 
            filter(! taxid %in% parent_taxid) # For each descendent, remove ancestors whose descendents are also ancestors
        taxids_extra_single <- taxids_extra_filtered %>% group_by(descendent_taxid) %>%
            filter(n() == 1) %>%
            select(taxid = descendent_taxid, parent_taxid = taxid) # If only one surviving ancestor, take that one
        taxids_extra_multi <- taxids_extra_filtered %>%
            filter(!descendent_taxid %in% taxids_extra_single$taxid) %>%
            mutate(rank_index = match(rank, ranks),
                   rank_index = replace_na(rank_index, 0)) %>%
            group_by(descendent_taxid) %>%
            filter(rank_index == min(rank_index)) %>%
            filter(taxid == max(taxid)) %>%
            filter(n() == 1) %>%
            select(taxid = descendent_taxid, parent_taxid = taxid) # If still multiple possible ancestors, pick the lowest rank
        taxids_extra_out <- bind_rows(taxids_extra_single, taxids_extra_multi)

        # Bind rows and write
        taxids_all_out <- bind_rows(hv_db_out, taxids_extra_out) %>%
            arrange(taxid) %>%
            mutate(rank = replace_na(rank, "no rank"),
                   name = replace_na(name, "unknown virus"))
        write_tsv(taxids_all_out, out_path)
        '''
