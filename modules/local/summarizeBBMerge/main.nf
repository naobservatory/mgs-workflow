process SUMMARIZE_BBMERGE {
    label "base"
    label "single"
    input:
        tuple val(sample), path(merged)
    output:
        path("${sample}_bbmerge_merged_summary.txt")
    script:
        """
        zcat ${merged} | awk '
            BEGIN {print "seq_id\tfrag_length"}
            NR % 4 == 1 {seq_id = \$1}
            NR % 4 == 2 {print seq_id "\t" length(\$0)}
        ' > ${sample}_bbmerge_merged_summary.txt
        echo Processing complete. Output saved to ${sample}_bbmerge_merged_summary.txt
        """
}