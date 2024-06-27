// Make taxid/name/rank DB for all viruses in NCBI taxonomy
process MAKE_TOTAL_VIRUS_DB {
    label "single"
    label "R"
    input:
        path(human_virus_db)
        path(taxonomy_nodes)
        path(taxonomy_names)
    output:
        path("total-virus-db.tsv.gz"), emit: db
    shell:
        '''
        #!/usr/bin/env Rscript

        # Setup
        library(tidyverse)
        nodes_path <- "!{taxonomy_nodes}"
        names_path <- "!{taxonomy_names}"
        hv_path    <- "!{human_virus_db}"
        out_path   <- "total-virus-db.tsv.gz"

        # Import data (NB: convoluted method for names file to avoid problems from unpaired strings)
        nodes <- read_tsv(nodes_path, col_names = FALSE) %>%
            transmute(taxid = X1, parent_taxid = X3, rank = X5)
        names <- names_path %>% read_file %>% gsub("\\"", "", .) %>%
            read_tsv(col_names = FALSE) %>%
            transmute(taxid = X1, name = X3, alt_name = X5, type = X7)
        hv    <- read_tsv(hv_path, col_names = TRUE)
        print(sum(is.na(hv)))

        # Get complete list of viral taxids
        taxids_in <- c(10239, hv$taxid) # 10239 = Viruses
        taxids_out <- taxids_in
        taxids_new <- filter(nodes, parent_taxid %in% taxids_out, !taxid %in% taxids_out) %>%
            pull(taxid) %>% sort
        while (length(taxids_new) > 0){
            taxids_out <- c(taxids_out, taxids_new) %>% sort
            taxids_new <- filter(nodes, parent_taxid %in% taxids_out, !taxid %in% taxids_out) %>%
                pull(taxid) %>% sort
        }
        print(sum(is.na(taxids_out)))

        # Generate nonhuman taxid DB from nodes and names
        nhv_db <- filter(nodes, taxid %in% taxids_out, ! taxid %in% hv$taxid) %>%
            left_join(names, by="taxid") %>% select(-alt_name)
        nhv_db_scinames <- filter(nhv_db, type == "scientific name") %>%
            group_by(taxid) %>% filter(row_number() == 1) %>% select(-type)
        nhv_db_other <- filter(nhv_db, ! taxid %in% nhv_db_scinames$taxid) %>%
            group_by(taxid) %>% filter(row_number() == 1) %>% select(-type)
        nhv_db_out <- bind_rows(nhv_db_scinames, nhv_db_other) %>%
            select(taxid, name, rank, parent_taxid) %>% arrange(taxid)
        print(sum(is.na(nhv_db)))
        print(sum(is.na(nhv_db_scinames)))
        print(sum(is.na(nhv_db_other)))
        print(sum(is.na(nhv_db_out)))

        # Bind rows and write
        virus_db_out <- bind_rows(nhv_db_out, hv) %>%
            arrange(taxid) %>%
            mutate(rank = replace_na(rank, "no rank"),
                   name = replace_na(name, "unknown virus"))
        print(sum(is.na(virus_db_out)))
        write_tsv(virus_db_out, out_path)
        '''
}
