// ABOUTME: Generic process to filter TSV files by column value, keeping rows that match specified criteria.
// ABOUTME: Takes column name and value as parameters, making it reusable across different filtering scenarios.
process FILTER_TSV_COLUMN_BY_VALUE {
    label "python"
    label "single"
    input:
        tuple val(sample), path(input_tsv)
        val(column_name)
        val(filter_value)
    output:
        tuple val(sample), path("${sample}_filtered.tsv.gz"), emit: output
        tuple val(sample), path("input_${input_tsv}"), emit: input
    script:
        """
        # Set strict error handling
        set -o pipefail
        out=${sample}_filtered.tsv.gz
        # Filter the TSV to keep only rows where column matches value
        zcat ${input_tsv} | filter_tsv_column_by_value.py -c "${column_name}" -v "${filter_value}" -o \${out}
        # Link input to output for testing
        ln -s ${input_tsv} input_${input_tsv}
        """
}