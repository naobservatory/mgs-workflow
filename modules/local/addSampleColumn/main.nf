// Add a column to a TSV with sample ID
process ADD_SAMPLE_COLUMN {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(tsv)
        val(sample_column)
        val(label)
    output:
        path("${sample}_${label}_out.tsv.gz"), emit: output
        path("${sample}_${label}_in.tsv.gz"), emit: input
    shell:
        '''
        add_sample_column.py !{tsv} !{sample} !{sample_column} !{sample}_!{label}_out.tsv.gz
        # Link input files for testing
        ln -s !{report} !{sample}_!{label}_in.tsv.gz
        '''
}
