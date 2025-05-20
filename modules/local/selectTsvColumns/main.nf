/*
Given a TSV column with an initial header line,
return a new TSV containing only the specified columns.
*/

process SELECT_TSV_COLUMNS {
    label "python"
    label "single"
    input:
        tuple val(sample), path(tsv) // Input TSV
        val(fields) // Comma-separated list of fields to select
    output:
        tuple val(sample), path("selected_${tsv}"), emit: output
        tuple val(sample), path("input_${tsv}"), emit: input
    shell:
        '''
        select_tsv_columns.py -i !{tsv} -o selected_!{tsv} -f !{fields}
        # Link input to output for testing
        ln -s !{tsv} input_!{tsv}
        '''
}
