// Process & filter BLAST output into TSV
process FILTER_BLAST {
    label "tidyverse"
    label "single_cpu_32GB_memory"
    input:
        path(blast_hits)
        val(max_rank)
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
        blast_hits_highrank <- blast_hits_ranked %>% filter(rank <= !{max_rank})
        # Write output
        write_tsv(blast_hits_highrank, "blast_hits_filtered.tsv.gz")
        '''
}

// Streamed version (final filter only, earlier filtering handled by FILTER_TSV)
process FILTER_BLAST_STREAMED {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(blast_hits_sorted) // Must be sorted on query ID (ascending) and bitscore (descending)
        val(max_rank) // Maximum bitscore rank to keep
        val(min_frac) // Minimum bitscore to retain (as a fraction of the best bitscore for the query)
    output:
        tuple val(sample), path("${sample}_blast_filtered.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_blast_in.tsv.gz"), emit: input
    shell:
        '''
        # Run script
        filter_blast.py -i !{blast_hits_sorted} -o !{sample}_blast_filtered.tsv.gz -r !{max_rank} -f !{min_frac}
        # Link input to output for testing
        ln -s !{blast_hits_sorted} !{sample}_blast_in.tsv.gz
        '''
}
