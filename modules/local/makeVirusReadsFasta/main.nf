// Extract FASTA from virus read TSV
process MAKE_VIRUS_READS_FASTA {
    label "tidyverse"
    label "single_cpu_16GB_memory"
    input:
        path(virus_hits_db)
    output:
        path("virus_hits_{1,2}.fasta.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        tab <- read_tsv("!{virus_hits_db}", col_names = TRUE, show_col_types = FALSE) %>%
            mutate(seq_head = paste0(">", seq_id),
                   header1 = paste0(seq_head, "_1"),
                   header2 = paste0(seq_head, "_2"))
        fasta_1_tab <- select(tab, header=header1, seq=query_seq_fwd)
        fasta_2_tab <- select(tab, header=header2, seq=query_seq_rev)
        fasta_1_out <- do.call(paste, c(fasta_1_tab, sep="\n")) %>% paste(collapse="\n")
        fasta_2_out <- do.call(paste, c(fasta_2_tab, sep="\n")) %>% paste(collapse="\n")
        out_file_1  <- gzfile("virus_hits_1.fasta.gz", "w")
        out_file_2  <- gzfile("virus_hits_2.fasta.gz", "w")
        write(fasta_1_out, out_file_1)
        write(fasta_2_out, out_file_2)
        close(out_file_1)
        close(out_file_2)
        '''
}
