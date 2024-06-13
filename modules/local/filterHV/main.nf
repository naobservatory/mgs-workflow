// Filter HV reads by alignment score
process FILTER_HV {
    label "tidyverse"
    label "single"
    input:
        path(hv_hits)
        val(score_threshold)
    output:
        path("hv_hits_putative_filtered.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        score_threshold <- !{score_threshold}
        data <- read_tsv("!{hv_hits}", col_names = TRUE, show_col_types = FALSE)
        filtered <- mutate(data, hit_hv = as.logical(!is.na(str_match(encoded_hits, paste0(" ", as.character(taxid), ":"))))) %>%
            mutate(adj_score_fwd = replace_na(adj_score_fwd, 0), adj_score_rev = replace_na(adj_score_rev, 0)) %>%
            filter((!classified) | assigned_hv) %>% 
            filter(adj_score_fwd > score_threshold | adj_score_rev > score_threshold | assigned_hv)
        print(dim(data))
        print(dim(filtered))
        write_tsv(filtered, "hv_hits_putative_filtered.tsv.gz")
        '''
}
