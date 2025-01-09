// Sort a gzipped TSV by a specified column header
// TODO: Expand to handle plaintext TSVs
process SORT_TSV {
    label "core_utils"
    label "single"
    input:
        tuple val(sample), path(tsv)
        val(sort_field)
        val(label)
    output:
        tuple val(sample), path("${sample}_${label}_sorted_${sort_field}.tsv.gz"), emit: sorted
        tuple val(sample), path("${sample}_${label}_in.tsv.gz"), emit: input
    shell:
        '''
        set -exuo pipefail
        # Extract header and determine column index to sort
        set +o pipefail
        HEADER=$(zcat !{tsv} | head -n 1)
        set -o pipefail
        # Identify index of sort field, and error out if it's absent
        IFS=$'\t' read -r -a header_fields <<< "$HEADER"
        COL_INDEX=0
        for i in "${!header_fields[@]}"; do
            if [[ "${header_fields[$i]}" == "!{sort_field}" ]]; then
                COL_INDEX=$(( i + 1 ))
                break
            fi
        done
        if [[ "${COL_INDEX}" -eq 0 ]]; then
            echo "Error: Could not find column '!{sort_field}' in the input header." >&2
            exit 1
        fi
        echo "Found '!{sort_field}' in column ${COL_INDEX}. Sorting by column ${COL_INDEX}."
        # Perform sorting while keeping header at top
        OUTPUT="!{sample}_!{label}_sorted_!{sort_field}.tsv.gz"
        set +o pipefail
        zcat !{tsv} | head -n 1 | gzip > ${OUTPUT}
        set -o pipefail
        zcat !{tsv} | tail -n +2 | sort -t $'\t' -k"${COL_INDEX}","${COL_INDEX}" | gzip >> ${OUTPUT}
        # Link input to output for testing
        ln -s !{tsv} !{sample}_!{label}_in.tsv.gz
        '''
}
