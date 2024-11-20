// Collapse filtered BLAST results across read pairs
process PAIR_BLAST {
    label "tidyverse"
    label "single"
    input:
        path(blast_filtered_fwd)
        path(blast_filtered_rev)
    output:
        path("blast_hits_paired.tsv.gz")
    shell:
        '''
        #!/usr/bin/env Rscript
        library(tidyverse)
        # Set input paths
        path_fwd <- "!{blast_filtered_fwd}"
        path_rev <- "!{blast_filtered_rev}"
        # Import data
        ctypes <- cols(qseqid="c", sseqid="c", sgi="c", staxid="i", qlen="i", evalue="d",
            bitscore="i", qcovs="i", length="i", pident="d", mismatch="i", gapopen="i",
            sstrand="c", qstart="i", qend="i", sstart="i", send="i", rank="i")
        hits_fwd <- read_tsv(path_fwd, col_types = ctypes)
        hits_rev <- read_tsv(path_rev, col_types = ctypes)
        # Join on query ID and subject taxid
        hits_pair <- full_join(hits_fwd, hits_rev, by=c("qseqid", "staxid"),
                suffix = c("_fwd", "_rev")) %>%
            mutate(n_reads = as.numeric(!is.na(sseqid_fwd))+as.numeric(!is.na(sseqid_rev)),
                bitscore_max = pmax(bitscore_fwd, bitscore_rev),
                bitscore_min = pmin(bitscore_fwd, bitscore_rev))
        # Write output
        write_tsv(hits_pair, "blast_hits_paired.tsv.gz")
        '''
}
