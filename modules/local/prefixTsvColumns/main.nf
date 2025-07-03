// Add prefix to TSV column headers with include/exclude mode
process PREFIX_TSV_COLUMNS {
    label "python"
    label "single"
    input:
        tuple val(sample), path(tsv)
        val(prefix)
        val(columns)    // comma-separated list of columns
        val(mode)       // "include" or "exclude"
    output:
        tuple val(sample), path("prefixed_${tsv}"), emit: output
        tuple val(sample), path("input_${tsv}"), emit: input
    shell:
        '''
        prefix_tsv_columns.py -i !{tsv} -o prefixed_!{tsv} -p "!{prefix}" -c "!{columns}" -m !{mode}
        # Link input to output for testing
        ln -s !{tsv} input_!{tsv}
        '''
}
