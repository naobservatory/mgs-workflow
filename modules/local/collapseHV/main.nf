// Collapse separate read pair entries in HV DB
process COLLAPSE_HV {
    label "tidyverse"
    cpus 1
    memory "16.GB"
    input:
        path(hv_hits_filtered)
        path(merged_bbmerge_paired_results)
        path(merged_bbmerge_concat_results)
        path(merged_alignment_dups)
    output:
        path("hv_hits_putative_collapsed.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        rmax <- function(x){
            if (all(is.na(x))) return(NA)
            return(max(x, na.rm = TRUE))
        }
        collapse <- function(x) ifelse(all(x == x[1]), x[1], paste(x, collapse="/"))
        hits_filtered <- read_tsv("!{hv_hits_filtered}", col_names = TRUE, show_col_types = FALSE)
        print(dim(hits_filtered))
        reads_collapsed <- hits_filtered %>% group_by(seq_id) %>% summarize(
            sample = collapse(sample), genome_id = collapse(genome_id),
            taxid_all = collapse(as.character(taxid)),
            taxid = taxid[1],
            best_alignment_score_fwd = rmax(best_alignment_score_fwd),
            best_alignment_score_rev = rmax(best_alignment_score_rev),
            query_len_fwd = rmax(query_len_fwd), query_seq_fwd = query_seq_fwd[!is.na(query_seq_fwd)][1],
            query_len_rev = rmax(query_len_rev), query_seq_rev = query_seq_rev[!is.na(query_seq_rev)][1],
            classified = rmax(classified), assigned_name = collapse(assigned_name),
            assigned_taxid_all = collapse(as.character(assigned_taxid)),
            assigned_taxid = assigned_taxid[1],
            assigned_hv = rmax(assigned_hv), hit_hv = rmax(hit_hv), encoded_hits = collapse(encoded_hits),
            adj_score_fwd = rmax(adj_score_fwd), adj_score_rev = rmax(adj_score_rev)
            ) %>% mutate(adj_score_max = pmax(adj_score_fwd, adj_score_rev))
        print(dim(reads_collapsed))
      
        bbmerge_merged <- read_tsv("!{merged_bbmerge_paired_results}", col_names = TRUE, show_col_types = FALSE) %>% 
            mutate(bbmerge_merge_status = 1) %>% rename(bbmerge_frag_length=frag_length)
        

        dedup <- read_tsv("!{merged_bbmerge_concat_results}", col_names = TRUE, show_col_types = FALSE) %>% 
            rename(merged_concatenated_sequence=sequence, clumpify_dupcount=counts)

        bbtools_summary <- left_join(dedup, bbmerge_merged, by = "seq_id") %>%
            mutate(bbmerge_merge_status = ifelse(is.na(bbmerge_merge_status), 0, bbmerge_merge_status),
            seq_id = str_remove(seq_id, "@")
            ) 
        
        alignment_dups <- read_tsv("!{merged_alignment_dups}", col_names = TRUE, show_col_types = FALSE) %>%
            rename(seq_id = query_name, 
            bowtie2_frag_length = fragment_length, 
            bowtie2_exemplar = exemplar,
            bowtie2_dupcount = dup_count)

        reads_collapsed_with_bbtools_summary <- left_join(reads_collapsed, bbtools_summary, by = "seq_id") %>%
            left_join(alignment_dups, by = "seq_id")

        write_tsv(reads_collapsed_with_bbtools_summary, "hv_hits_putative_collapsed.tsv.gz")
        '''
}
