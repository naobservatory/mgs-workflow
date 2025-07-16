// Creates a new column with conditional values based on checking another column
process CREATE_CONDITIONAL_COLUMN {
    label "single"
    label "coreutils_gzip_gawk"
    input:
        tuple val(sample), path(tsv)
        val(chk_col)
        val(match_val)
        val(if_col)
        val(else_col)
        val(new_hdr)
    output:
        tuple val(sample), path("adjusted_${new_hdr}_${tsv}"), emit: tsv
    script:
        def isGzipped = tsv.toString().endsWith(".gz")
        def extractCmd = isGzipped ? "zcat" : "cat"
        def outputFile = "adjusted_${new_hdr}_${tsv}"
        def compressCmd = isGzipped ? "| gzip" : ""
        """
        set -eou pipefail
        ${extractCmd} ${tsv} | awk -v chk_col="${chk_col}" \\
            -v match_val="${match_val}" \\
            -v if_col="${if_col}" \\
            -v else_col="${else_col}" \\
            -v new_hdr="${new_hdr}" \\
            'BEGIN {
                FS = OFS = "\\t"
                error_occurred = 0
            }
            # Read the header line, check if the columns exist, and add the new column
            NR == 1 {
                if (NF > 0) {
                    for (i = 1; i <= NF; i++) {
                        if (\$i == chk_col)   chk_idx   = i
                        if (\$i == if_col)    if_idx    = i
                        if (\$i == else_col)  else_idx  = i
                    }
                    if (!(chk_idx && if_idx && else_idx)) {
                        printf("ERROR: could not find all requested columns in header\\n chk_col=%s (idx=%d) if_col=%s (idx=%d) else_col=%s (idx=%d)\\n", chk_col, chk_idx, if_col, if_idx, else_col, else_idx) > "/dev/stderr"
                        error_occurred = 1
                        exit 1
                    }
                    print \$0, new_hdr
                }
                next
            }
            # Read the rest of the file, and add the new column
            {
                \$(NF + 1) = ( \$chk_idx == match_val ? \$(if_idx) : \$(else_idx) )
                print
            }
            END {
                # Make sure to exit with a non-zero status if a header exists, but the columns could not be found
                if (error_occurred) {
                    exit 1
                }
                if (NR == 0 || NR == 1) {
                    exit 0
                }
            }' ${compressCmd} > ${outputFile}
        """
}
