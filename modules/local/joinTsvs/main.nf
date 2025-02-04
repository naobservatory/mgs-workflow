// Join two sorted TSVs on a specified column header
process JOIN_TSVS {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(tsv1), path(tsv2)
        val(join_field)
        val(join_type) // inner, left, right, outer, strict
        val(label)
    output:
        tuple val(sample), path("${sample}_${label}_${join_type}_joined_${join_field}.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_${label}_in_{1,2}.tsv.gz"), emit: input
    shell:
        '''
        out=!{sample}_!{label}_!{join_type}_joined_!{join_field}.tsv.gz
        join_tsvs.py !{tsv1} !{tsv2} !{join_field} !{join_type} ${out}
        # Link input files to output for testing
        ln -s !{tsv1} !{sample}_!{label}_in_1.tsv.gz
        ln -s !{tsv2} !{sample}_!{label}_in_2.tsv.gz
        '''
}
