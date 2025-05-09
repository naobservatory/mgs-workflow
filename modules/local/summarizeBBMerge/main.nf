process SUMMARIZE_BBMERGE {
    label "coreutils_gzip_gawk"
    label "single"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_bbmerge_summary.tsv.gz"), emit: summary
        tuple val(sample), path("${sample}_bbmerge_summary_in_merged.fastq.gz"), emit: input
    shell:
        '''
        # Check if input file is empty
        if [[ ! -s !{reads[0]} ]] || [[ $(zcat !{reads[0]} | head -c1 | wc -c) -eq 0 ]]; then
            echo "Warning: Input file is empty. Creating empty summary with header only."
            echo -e "seq_id\tbbmerge_frag_length" | gzip -c > !{sample}_bbmerge_summary.tsv.gz
        else
            # Process non-empty file
            zcat !{reads[0]} | awk '
                BEGIN {print "seq_id\tbbmerge_frag_length"}
                NR % 4 == 1 {
                    sub(/^@/, "", $1)
                    seq_id = $1
                }
                NR % 4 == 2 {print seq_id "\t" length($0)}
            ' | gzip -c > !{sample}_bbmerge_summary.tsv.gz
        fi
        
        echo "Processing complete. Output saved to !{sample}_bbmerge_summary.tsv.gz"
        ln -s !{reads[0]} !{sample}_bbmerge_summary_in_merged.fastq.gz
        '''
}
