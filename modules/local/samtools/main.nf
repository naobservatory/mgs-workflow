// Return reads that did not align to reference as FASTQ (streamed version)
process SAMTOOLS_FILTER {
    label "samtools"
    label "small"
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

// Return aligned and unaligned reads separately as FASTQs
process SAMTOOLS_SEPARATE {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_samtools_fail.fastq.gz"), emit: nomatch
        tuple val(sample), path("${sample}_${suffix}_samtools_pass.fastq.gz"), emit: match
        tuple val(sample), path("${sample}_in.sam"), emit: input
    shell:
        '''
        # Define output
        of=!{sample}_!{suffix}_samtools_fail.fastq.gz
        op=!{sample}_!{suffix}_samtools_pass.fastq.gz
        of_var="fastq -n -F 4"
        op_var="fastq -n -f 4"
        # Execute samtools
        samtools ${of_var} ${in} | gzip > ${of}
        samtools ${op_var} ${in} | gzip > ${op}
        # Link input to output for testing
        ln -s !{sam} !{sample}_in.sam
        '''
}