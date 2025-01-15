// Rename a field in a gzipped TSV header
// TODO: Expand to handle plaintext TSVs
// TODO: Expand to rename multiple header fields at once
process REHEAD_TSV {
    label "core_utils"
    label "single"
    input:
        tuple val(sample), path(tsv)
        val(old_field)
        val(new_field)
        val(label)
    output:
        tuple val(sample), path("${sample}_${label}_rename_${old_field}_${new_field}.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_${label}_in.tsv.gz"), emit: input
    shell:
        '''
        set -exuo pipefail
        gunzip -c !{tsv} | {
            # Read header line
            read -r HEADER
            # Parse into a tab-separated array
            IFS=$'\t' read -r -a fields <<< "${HEADER}"
            # Replace old field with new field
            FOUND=0
            for i in "${!fields[@]}"; do
                if [[ "${fields[$i]}" == "!{old_field}" ]]; then
                    fields[$i]="!{new_field}"
                    FOUND=1
                fi
            done
            # Error out if we never found old_field
            if [[ "${FOUND}" -eq 0 ]]; then
                echo "ERROR: Target field not found in header: '!{old_field}'" >&2
                exit 1
            fi
            # Print the modified header
            printf "%s\n" "$(IFS=$'\t'; echo "${fields[*]}")"
            # Copy the rest of the lines verbatim
            cat
        } | gzip -c > !{sample}_!{label}_rename_!{old_field}_!{new_field}.tsv.gz
        # Link input to output for testing
        ln -s !{tsv} !{sample}_!{label}_in.tsv.gz
        '''
}
