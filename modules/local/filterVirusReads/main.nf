// Filter virus reads by alignment score and assignment status
process FILTER_VIRUS_READS {
    label "tidyverse"
    cpus 1
    memory "16.GB"
    input:
        path(virus_hits)
        val(score_threshold)
    output:
        path("virus_hits_putative_filtered.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        score_threshold <- !{score_threshold}
        data <- read_tsv("!{virus_hits}", col_names = TRUE, show_col_types = FALSE)
        filtered <- mutate(data, kraken_hit_host_virus = as.logical(!is.na(str_match(kraken_encoded_hits, paste0(" ", as.character(taxid), ":"))))) %>%
            mutate(adj_score_fwd = replace_na(adj_score_fwd, 0), adj_score_rev = replace_na(adj_score_rev, 0)) %>%
            filter((!kraken_classified) | kraken_assigned_host_virus > 0) %>% # Remove reads that are assigned to a non-host-virus taxon
            filter(adj_score_fwd > score_threshold | adj_score_rev > score_threshold | kraken_assigned_host_virus == 1)
        print(dim(data))
        print(dim(filtered))
        write_tsv(filtered, "virus_hits_putative_filtered.tsv.gz")
        '''
}
