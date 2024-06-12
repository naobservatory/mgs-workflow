// Combine MultiQC summary files across workflow stages
process MERGE_MULTIQC {
    label "tidyverse"
    label "single"
    input:
        path(basic_stats_tsvs)
        path(adapter_tsvs)
        path(base_quality_tsvs)
        path(sequence_quality_tsvs)
    output:
        path("qc_basic_stats.tsv.gz")
        path("qc_adapter_stats.tsv.gz")
        path("qc_quality_base_stats.tsv.gz")
        path("qc_quality_sequence_stats.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        # Import data
        in_paths_basic <- str_split("!{basic_stats_tsvs}", " ")[[1]]
        in_paths_adapt <- str_split("!{adapter_tsvs}", " ")[[1]]
        in_paths_qbase <- str_split("!{base_quality_tsvs}", " ")[[1]]
        in_paths_qseqs <- str_split("!{sequence_quality_tsvs}", " ")[[1]]
        tabs_basic <- lapply(in_paths_basic, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        tabs_adapt <- lapply(in_paths_adapt, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        tabs_qbase <- lapply(in_paths_qbase, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        tabs_qseqs <- lapply(in_paths_qseqs, function(t) read_tsv(t, col_names = TRUE, cols(.default="c")))
        # Bind rows
        tab_basic_out <- bind_rows(tabs_basic)
        tab_adapt_out <- bind_rows(tabs_adapt)
        tab_qbase_out <- bind_rows(tabs_qbase)
        tab_qseqs_out <- bind_rows(tabs_qseqs)
        # Write output
        write_tsv(tab_basic_out, "qc_basic_stats.tsv.gz")
        write_tsv(tab_adapt_out, "qc_adapter_stats.tsv.gz")
        write_tsv(tab_qbase_out, "qc_quality_base_stats.tsv.gz")
        write_tsv(tab_qseqs_out, "qc_quality_sequence_stats.tsv.gz")
        '''
}
