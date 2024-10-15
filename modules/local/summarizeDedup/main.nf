process SUMMARIZE_DEDUP {
    label "base"
    label "single"
    input:
        tuple val(sample), path(merged)
    output:
        path("${sample}_dedup_summary.txt")
    script:
        """
        zcat ${merged} | awk '
            BEGIN {print "seq_id\tmerged_concatenated_sequence\tclumpify_dupcount"}
            NR % 4 == 1 {
                seq_id = \$1
                counts = 1
                if (match(\$0, /copies=([0-9]+)/)) {
                    counts = substr(\$0, RSTART+7, RLENGTH-7)
                }
            }
            NR % 4 == 2 {
                print seq_id "\t" \$0 "\t" counts
            }
        ' > ${sample}_dedup_summary.txt
        echo Processing complete. Output saved to ${sample}_dedup_summary.txt
        """
}