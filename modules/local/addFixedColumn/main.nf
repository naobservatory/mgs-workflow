// Add a column to a TSV with a specified name and value
process ADD_FIXED_COLUMN {
    label "biopython"
    label "single"
    input:
        path(tsv)
        val(column)
        val(value)
        val(label)
    output:
        path("${label}_labeled.tsv.gz"), emit: output
        path("${label}_in.tsv.gz"), emit: input
    shell:
        '''
        add_fixed_column.py !{tsv} !{column} !{value} !{label}_labeled.tsv.gz
        # Link input files for testing
        ln -s !{tsv} !{label}_in.tsv.gz
        '''
}
