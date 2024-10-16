// Collapse separate read pair entries in HV DB
process ADD_FRAG_DUP_TO_HV {
    label "tidyverse"
    cpus 1
    memory "16.GB"
    input:
        path(collapsed_ch)
        path(merged_bbmerge_results)
        path(merged_dedup_results)
    output:
        path("hv_hits_putative_collapsed.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)

        reads_collapsed <- read_tsv("!{collapsed_ch}", col_names = TRUE, show_col_types = FALSE)
      
        bbmerge <- read_tsv("!{merged_bbmerge_results}", col_names = TRUE, show_col_types = FALSE) %>% 
            mutate(bbmerge_merge_status = 1)

        dedup <- read_tsv("!{merged_dedup_results}", col_names = TRUE, show_col_types = FALSE)

        # Print column headers for bbmerge
        cat("BBMerge column headers:\n")
        print(colnames(bbmerge))

        # Print column headers for dedup
        cat("\nDedup column headers:\n")
        print(colnames(dedup))

        bbtools_summary <- left_join(dedup, bbmerge, by = "seq_id") %>%
            mutate(
                bbmerge_merge_status = ifelse(is.na(bbmerge_merge_status), 0, bbmerge_merge_status),
                seq_id = str_remove(seq_id, "@")
            ) 
        
        alignment_dup <- read_tsv("!{merged_alignment_dups}", col_names = TRUE, show_col_types = FALSE) %>%
            rename(seq_id = query_name, 
            bowtie2_frag_length = fragment_length, 
            bowtie2_exemplar = exemplar,
            bowtie2_dupcount = dup_count)
        
        reads_collapsed_with_bbtools_summary <- left_join(reads_collapsed, bbtools_summary, by = "seq_id") %>%
            left_join(alignment_dup, by = "seq_id")

        write_tsv(reads_collapsed_with_bbtools_summary, "hv_hits_putative_collapsed.tsv.gz")
        '''
}
