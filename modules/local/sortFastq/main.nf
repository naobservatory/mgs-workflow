// Sort a gzipped FASTQ file based on header sequences
process SORT_FASTQ {
    label "coreutils"
    label "single"
    input:
        tuple val(sample), path(fastq) // Interleaved or single-end
    output:
        tuple val(sample), path("sorted_${fastq}"), emit: output
        tuple val(sample), path("input_${fastq}"), emit: input
    shell:
        '''
        set -euo pipefail
        zcat !{fastq} | paste - - - - | sort -k1,1 | \\
            tr '\\t' '\\n' | gzip > sorted_!{fastq}
        # Link input to output for testing
        ln -s !{fastq} input_!{fastq}
        '''
}
