// 8.3. Summarize taxonomic composition from Bracken and MultiQC output
process SUMMARIZE_COMPOSITION {
    label "tidyverse"
    label "single"
    input:
        path(bracken)
        path(basic_stats)
        val(stage)
        val(classified_subset)
    output:
        path("taxonomic_composition.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        # Import data
        bracken <- read_tsv("!{bracken}", show_col_types = FALSE)
        basic   <- read_tsv("!{basic_stats}", show_col_types = FALSE)
        # Compute reads assigned to specific taxa
        read_comp_taxa_raw <- bracken %>% select(name, sample, new_est_reads) %>%
          mutate(name = tolower(name)) %>%
          pivot_wider(id_cols = "sample", names_from = "name", values_from = "new_est_reads",
                      values_fill = 0)
        read_comp_taxa <- select(read_comp_taxa_raw, sample)
        if ("bacteria" %in% colnames(read_comp_taxa_raw)) read_comp_taxa[["n_bacterial"]] <- read_comp_taxa_raw[["bacteria"]]
        if ("archaea" %in% colnames(read_comp_taxa_raw)) read_comp_taxa[["n_archaeal"]] <- read_comp_taxa_raw[["archaea"]]
        if ("viruses" %in% colnames(read_comp_taxa_raw)) read_comp_taxa[["n_viral"]] <- read_comp_taxa_raw[["viruses"]]
        if ("eukaryota" %in% colnames(read_comp_taxa_raw)) read_comp_taxa[["n_human"]] <- read_comp_taxa_raw[["eukaryota"]]
        # Compute reads surviving at each included preprocessing stage
        total_assigned <- bracken %>% group_by(sample) %>% summarize(
          name = "Total",
          kraken_assigned_reads = sum(kraken_assigned_reads),
          added_reads = sum(added_reads),
          new_est_reads = sum(new_est_reads),
          fraction_total_reads = sum(fraction_total_reads)
        )
        read_counts_preproc <- basic %>% select(sample, stage, n_read_pairs) %>%
          mutate(n_read_pairs = n_read_pairs * !{classified_subset}) %>%
          pivot_wider(id_cols = c("sample"), names_from="stage", values_from="n_read_pairs",
                      values_fill = 0) %>%
          inner_join(total_assigned %>% select(sample, new_est_reads), by = "sample") %>%
          rename(assigned = new_est_reads)
        if ("!{stage}" == "raw_concat"){
            read_comp_preproc <- transmute(read_counts_preproc, sample=sample, n_unassigned=raw_concat-assigned)
        } else if ("!{stage}" == "cleaned"){
            read_comp_preproc <- transmute(read_counts_preproc, sample=sample, n_filtered=raw_concat-cleaned,
                n_unassigned = cleaned-assigned)
        } else if ("!{stage}" == "dedup"){
            read_comp_preproc <- transmute(read_counts_preproc, sample=sample, n_filtered=raw_concat-cleaned,
                n_duplicate=cleaned-dedup, n_unassigned=dedup-assigned)
        } else if ("!{stage}" == "ribo_initial"){
            read_comp_preproc <- transmute(read_counts_preproc, sample=sample, n_filtered=raw_concat-cleaned,
                n_duplicate=cleaned-dedup, n_ribosomal=dedup-ribo_initial, n_unassigned=ribo_initial-assigned)
        } else if ("!{stage}" == "ribo_secondary"){
            read_comp_preproc <- transmute(read_counts_preproc, sample=sample, n_filtered=raw_concat-cleaned,
                n_duplicate=cleaned-dedup, n_ribosomal=dedup-ribo_secondary, n_unassigned=ribo_secondary-assigned)
        }
        # Compute overall composition
        read_comp <- full_join(read_comp_preproc, read_comp_taxa, by="sample")
        read_comp_long <- pivot_longer(read_comp, -(sample), names_to = "classification",
                                       names_prefix = "n_", values_to = "n_reads") %>%
          mutate(classification = fct_inorder(str_to_sentence(classification))) %>%
          group_by(sample) %>% mutate(p_reads = n_reads/sum(n_reads))
        # Write output
        write_tsv(read_comp_long, "taxonomic_composition.tsv.gz")
        '''
}
