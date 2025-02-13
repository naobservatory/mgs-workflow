process MARK_ALIGNMENT_DUPLICATES {
    label "single"
    label "coreutils"
    input:
        tuple val(sample), path(tsv)
        val(fuzzy_match)
    output:
        tuple val(sample), path("${sample}_duplicate_reads.tsv.gz"), emit: output
        tuple val(sample), path("input_${tsv}"), emit: input
    shell:
    '''
    mark_duplicates "!{tsv}" "!{sample}_duplicate_reads.tsv.gz" !{fuzzy_match}
    ln -s !{tsv} input_!{tsv} # Link output to input for testing
    '''
}
