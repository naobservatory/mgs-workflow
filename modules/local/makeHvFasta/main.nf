// Extract FASTA from HV sequence TSV
process MAKE_HV_FASTA {
    label "tidyverse"
    cpus 1
    memory "16.GB"
    input:
        path(hv_hits_collapsed)
    output:
        path("hv_hits_putative_{1,2}.fasta.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        tab <- read_tsv("!{hv_hits_collapsed}", col_names = TRUE, show_col_types = FALSE) %>%
            mutate(seq_head = paste0(">", seq_id),
                   header1 = paste0(seq_head, "_1"),
                   header2 = paste0(seq_head, "_2"))
        fasta_1_tab <- select(tab, header=header1, seq=query_seq_fwd)
        fasta_2_tab <- select(tab, header=header2, seq=query_seq_rev)
        fasta_1_out <- do.call(paste, c(fasta_1_tab, sep="\n")) %>% paste(collapse="\n")
        fasta_2_out <- do.call(paste, c(fasta_2_tab, sep="\n")) %>% paste(collapse="\n")
        out_file_1  <- gzfile("hv_hits_putative_1.fasta.gz", "w")
        out_file_2  <- gzfile("hv_hits_putative_2.fasta.gz", "w")
        write(fasta_1_out, out_file_1)
        write(fasta_2_out, out_file_2)
        close(out_file_1)
        close(out_file_2)
        '''
}
