// Extract sequence IDs from a FASTA or FASTQ file
process EXTRACT_FASTN_IDS {
    label "seqkit"
    label "single"
    label "testing" // Currently only used for testing
    input:
        tuple val(sample), path(fastn)
    output:
        tuple val(sample), path("${sample}_ids.txt"), emit: output
        tuple val(sample), path("input_${fastn}"), emit: input
    shell:
        '''
        seqkit seq -ni !{fastn} | sort -u > !{sample}_ids.txt
        ln -s !{fastn} input_!{fastn}
        '''
}
