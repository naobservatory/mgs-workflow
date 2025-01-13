process COMBINE_READ_COUNTS {
    label "coreutils_gzip_gawk"
    label "single"
    input:
        path(read_counts)
    output:
        path("read_counts.tsv.gz")
    shell:
        '''
        # Add header
        echo "sample,total_read_count" > read_counts.tsv
        # Add read counts
        cat !{read_counts.join(" ")} >> read_counts.tsv
        # Compress
        gzip read_counts.tsv
        '''
}