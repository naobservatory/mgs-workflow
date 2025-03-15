// Given a gzipped TSV of viral hits, extract the seq_ids and use them to subseq an interleaved FASTQ file
process EXTRACT_VIRAL_HITS_TO_FASTQ {
    label "seqtk"
    label "single_cpu_16GB_memory"
    input:
        path(hits_tsv)
        path(fastq)
    output:
        path("virus_hits_filtered.fastq.gz"), emit: fastq
        path("virus_hits_ids.txt.gz"), emit: ids
        path("virus_hits_in.tsv.gz"), emit: input_tsv
        path("virus_reads_in.fastq.gz"), emit: input_fastq
    shell:
        '''
        set -euo pipefail
        zcat !{hits_tsv} \\
            | awk -F'\t' '
                BEGIN { colIndex = -1 }
                NR == 1 {
                    for (i = 1; i <= NF; i++) { if ($i == "seq_id") { colIndex = i; break } }
                    if (colIndex == -1) {
                        print "ERROR: No column named 'seq_id' in header" > "/dev/stderr"
                        exit 1
                    } else {
                        print "'seq_id' column index: $colIndex" > "/dev/stderr"
                    }
                    next # Skip printing the header itself
                }
                { print $colIndex } # Print the seq_id column only from subsequent lines
            ' \\
            | tee >(gzip -c > virus_hits_ids.txt.gz) \\
            | seqtk subseq !{fastq} - \\
            | gzip -c > virus_hits_filtered.fastq.gz
        # Link input to output for testing
        ln -s !{hits_tsv} virus_hits_in.tsv.gz
        ln -s !{fastq} virus_reads_in.fastq.gz
        '''
}
