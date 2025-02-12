// Return reads that did not align to reference as FASTQ (streamed version)
process SAMTOOLS_FILTER {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}.fastq.gz"), emit: reads
        tuple val(sample), path("${sample}_in.sam"), emit: input
    shell:
        '''
        # Define output
        out=!{sample}_!{suffix}.fastq.gz
        var="fastq -n -f 4"
        # Execute samtools
        samtools ${var} !{sam} | gzip > ${out}
        # Link input to output for testing
        ln -s !{sam} !{sample}_in.sam
        '''
}