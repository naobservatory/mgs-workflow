// Process to add prefix to TSV column headers with include/exclude mode for flexible column selection
process PREFIX_TSV_COLUMNS {
    label "python"
    label "single"
    input:
        tuple val(sample), path(tsv)
        val(prefix)
        val(columns)    // comma-separated list of columns
        val(mode)       // "include" or "exclude"
    output:
        tuple val(sample), path("prefix_${mode}_${tsv}"), emit: output
        tuple val(sample), path("input_${mode}_${tsv}"), emit: input
    script:
        """
        # Set strict error handling
        set -o pipefail
        
        # Add prefix to columns based on mode
        prefix_tsv_columns.py -i ${tsv} -o prefix_${mode}_${tsv} -p "${prefix}" -c "${columns}" -m ${mode}
        
        # Link input to output for testing
        ln -s ${tsv} input_${mode}_${tsv}
        """
}
