// Rename a field in a gzipped TSV header
// TODO: Expand to handle plaintext TSVs
// TODO: Expand to rename multiple header fields at once
process REHEAD_TSV {
    label "core_utils"
    label "single"
    input:
        tuple val(sample), path(tsv)
        val(old_fields)
        val(new_fields)
        val(label)
    output:
        tuple val(sample), path("${sample}_${label}_renamed.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_${label}_in.tsv.gz"), emit: input
    shell:
        '''
        rehead_tsv.py !{tsv} !{old_fields} !{new_fields} !{sample}_!{label}_renamed.tsv.gz
        # Link input to output for testing
        ln -s !{tsv} !{sample}_!{label}_in.tsv.gz
        '''
}
