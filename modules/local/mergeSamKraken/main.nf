// Merge processed SAM and Kraken TSVs and compute length-normalized alignment scores
process MERGE_SAM_KRAKEN {
    label "tidyverse"
    label "single_cpu_16GB_memory"
    input:
        tuple val(sample), path(kraken_processed), path(sam_processed)
    output:
        path("${sample}_virus_hits_putative.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        sam <- read_tsv("!{sam_processed}", show_col_types = FALSE)
        krk <- read_tsv("!{kraken_processed}", show_col_types = FALSE)
        mrg <- sam %>% rename(seq_id = query_name) %>% inner_join(krk, by="seq_id")
        cat(nrow(sam), nrow(krk), nrow(mrg), "\n")
        cat(ncol(sam), ncol(krk), ncol(mrg), "\n")
        cat(names(mrg), "\n")
        cat(head(mrg$best_alignment_score_fwd))
        cat(head(mrg$query_len_fwd))
        cat(class(mrg$query_len_fwd))
        if (nrow(mrg)>0) {
            mrg <- mrg %>%
                mutate(adj_score_fwd = best_alignment_score_fwd/log(query_len_fwd),
                       adj_score_rev = best_alignment_score_rev/log(query_len_rev),
                       sample="!{sample}")
        }
        write_tsv(mrg, "!{sample}_virus_hits_putative.tsv.gz")
        '''
}
