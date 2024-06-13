// Collapse filtered BLAST results across read pairs
process PAIR_BLAST {
    label "tidyverse"
    label "single"
    input:
        path(blast_hits_filtered)
    output:
        path("blast_hits_paired.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        blast_hits_filtered_path <- "!{blast_hits_filtered}"
        blast_hits_filtered <- read_tsv(blast_hits_filtered_path, show_col_types = FALSE)
        # Summarize by read pair and taxid
        blast_hits_paired <- blast_hits_filtered %>%
            mutate(bitscore = as.numeric(bitscore)) %>%
            group_by(seq_id, staxid) %>%
            summarize(bitscore_max = max(bitscore), bitscore_min = min(bitscore),
                      n_reads = n(), .groups = "drop")
        # Write output
        write_tsv(blast_hits_paired, "hv_hits_blast_paired.tsv.gz")
        '''
}
