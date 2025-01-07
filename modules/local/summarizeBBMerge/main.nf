process SUMMARIZE_BBMERGE {
    label "coreutils_gzip_gawk"
    label "single"
    input:
        tuple val(sample), path(merged)
    output:
        path("${sample}_bbmerge_summary.txt")
    script:
        """
        zcat ${merged} | awk '
            BEGIN {print "seq_id\tbbmerge_frag_length"}
            NR % 4 == 1 {seq_id = \$1}
            NR % 4 == 2 {print seq_id "\t" length(\$0)}
        ' > ${sample}_bbmerge_summary.txt
        echo Processing complete. Output saved to ${sample}_bbmerge_summary.txt
        """
}

process SUMMARIZE_BBMERGE_STREAMED {
    label "coreutils_gzip_gawk"
    label "single"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_bbmerge_summary.tsv.gz"), emit: summary
        tuple val(sample), path("${sample}_bbmerge_summary_in_merged.fastq.gz"), emit: input
    script:
        """
        zcat ${reads[0]} | awk -v sample=${sample} '
            BEGIN {print "sample\tseq_id\tbbmerge_frag_length"}
            NR % 4 == 1 {seq_id = \$1}
            NR % 4 == 2 {print sample "\t" seq_id "\t" length(\$0)}
        ' | gzip -c > ${sample}_bbmerge_summary.tsv.gz
        echo Processing complete. Output saved to ${sample}_bbmerge_summary.txt
        ln -s ${reads[0]} ${sample}_bbmerge_summary_in_merged.fastq.gz
        """
}
