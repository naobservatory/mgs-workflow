// Add a column to a TSV with a specified name and value
process ADD_FIXED_COLUMN {
    label "python"
    label "single"
    input:
        path(tsv)
        val(column)
        val(value)
        val(label)
    output:
        path("labeled_${label}_${tsv}"), emit: output
        path("input_${tsv}"), emit: input
    shell:
        '''
        add_fixed_column.py !{tsv} !{column} !{value} labeled_!{label}_!{tsv}
        # Link input files for testing
        ln -s !{tsv} input_!{tsv}
        '''
}
