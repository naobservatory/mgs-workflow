// Mask low complexity FASTQ read regions. Only works on gzipped FASTQ files.
process MASK_FASTQ_READS {
    label "large"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
        val(window_size)
	    val(entropy)
    output:
        tuple val(sample), path("${sample}_masked.fastq.gz"), emit: masked
        tuple val(sample), path("${sample}_in.fastq.gz"), emit: input
    shell:
        '''
        set -eoux pipefail

        # Define input/output
        out=!{sample}_masked.fastq.gz

        # Define parameters
        par="window=!{window_size} entropy=!{entropy}"

        # If input is empty, create empty gzipped output (bbmask errors on empty input)
        if [[ -z $(zcat "${reads}" | head) ]]; then
            echo "Input file is empty. Creating empty output."
            echo -n | gzip > ${out}
        else
            # Execute with streaming approach
            echo "non-empty input"
            zcat -f !{reads} | bbmask.sh in=stdin.fastq out=stdout.fastq ${par} | gzip > ${out}
        
            # Check for empty output file without empty input
            if [[ -z ($(zcat "${out}" | head) ]]; then
                echo "Error: Output file is empty."
                exit 1
            fi
        fi

        # Link input to output for testing
        ln -s !{reads} !{sample}_in.fastq.gz
        '''
}