// Sort a gzipped TSV by a specified column header
// TODO: Expand to handle plaintext TSVs
process SORT_TSV {
    label "coreutils_file"
    label "single"
    input:
        tuple val(sample), path(tsv)
        val(sort_field)
    output:
        tuple val(sample), path("sorted_${sort_field}_${tsv}"), emit: sorted
        tuple val(sample), path("input_${tsv}"), emit: input
    shell:
        '''
        set -exuo pipefail
        # Handle gzipped input
        in_status=$(file -b $(readlink -f !{tsv}))
        echo "Input file status: ${in_status}"
        gzip_in=$(echo "${in_status}" | grep -c "gzip" || true)
        echo "Gzip count in file status: ${gzip_in}"
        if [[ ${gzip_in} -gt 0 ]]; then
            echo "Gzipped input - running zcat."
            fn=zcat
            fn2=gzip
        else
            echo "Non-gzipped input - running cat."
            fn=cat
            fn2=cat
        fi
        # Extract header and determine column index to sort
        set +o pipefail
        HEADER=$(${fn} !{tsv} | head -n 1)
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
            echo "ERROR: Could not find sort field in input header: '!{sort_field}'." >&2
            exit 1
        fi
        echo "Found '!{sort_field}' in column ${COL_INDEX}. Sorting by column ${COL_INDEX}."
        # Perform sorting while keeping header at top
        OUTPUT="sorted_!{sort_field}_!{tsv}"
        set +o pipefail
        ${fn} !{tsv} | head -n 1 | ${fn2} > ${OUTPUT}
        set -o pipefail
        ${fn} !{tsv} | tail -n +2 | sort -t $'\t' -k"${COL_INDEX}","${COL_INDEX}" | ${fn2} >> ${OUTPUT}
        # Link input to output for testing
        ln -s !{tsv} input_!{tsv}
        '''
}
