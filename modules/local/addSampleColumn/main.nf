// Add a column to a TSV with sample ID
process ADD_SAMPLE_COLUMN {
    label "python"
    label "single"
    input:
        tuple val(sample), path(tsv)
        val(sample_column)
        val(label)
    output:
        tuple val(sample), path("${sample}_${label}_labeled.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_${label}_in.tsv.gz"), emit: input
    shell:
        '''
        add_sample_column.py !{tsv} !{sample} !{sample_column} !{sample}_!{label}_labeled.tsv.gz
        # Link input files for testing
        ln -s !{tsv} !{sample}_!{label}_in.tsv.gz
        '''
}
