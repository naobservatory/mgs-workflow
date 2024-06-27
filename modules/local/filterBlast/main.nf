// Process & filter BLAST output into TSV
process FILTER_BLAST {
    label "tidyverse"
    cpus 1
    memory "${params.mem}"
    input:
        path(blast_hits)
    output:
        path("blast_hits_filtered.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        blast_hits_path <- "!{blast_hits}"
        blast_cols <- c("qseqid", "sseqid", "sgi", "staxid", "qlen", "evalue",
            "bitscore", "qcovs", "length", "pident", "mismatch", "gapopen",
            "sstrand", "qstart", "qend", "sstart", "send")
        blast_hits <- read_tsv(blast_hits_path, show_col_types = FALSE,
            col_names = blast_cols, col_types = cols(.default="c"))
        # Filter for best hit for each query/subject combination
        blast_hits_best <- blast_hits %>% group_by(qseqid, staxid) %>% 
          filter(bitscore == max(bitscore)) %>%
          filter(length == max(length)) %>% filter(row_number() == 1)
        # Rank hits for each query and filter for high-ranking hits
        blast_hits_ranked <- blast_hits_best %>% 
          group_by(qseqid) %>% mutate(rank = dense_rank(desc(bitscore)))
        blast_hits_highrank <- blast_hits_ranked %>% filter(rank <= 5) %>%
            mutate(read_pair = str_split(qseqid, "_") %>% sapply(nth, n=-1),
                   seq_id = str_split(qseqid, "_") %>% sapply(nth, n=1))
        # Write output
        write_tsv(blast_hits_highrank, "blast_hits_filtered.tsv.gz")
        '''
}
