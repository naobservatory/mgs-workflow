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
        cat !{sam} | samtools ${var} - | gzip > ${out}
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
        tuple val(sample), path("${sample}_${suffix}_samtools_match.fastq.gz"), emit: match
        tuple val(sample), path("${sample}_${suffix}_samtools_nomatch.fastq.gz"), emit: nomatch
        tuple val(sample), path("${sample}_in.sam"), emit: input
    shell:
        '''
        # Define output
        om=!{sample}_!{suffix}_samtools_match.fastq.gz
        on=!{sample}_!{suffix}_samtools_nomatch.fastq.gz
        om_var="fastq -n -F 4"
        on_var="fastq -n -f 4"
        # Execute samtools
        cat !{sam} | samtools ${om_var} - | gzip > ${om}
        cat !{sam} | samtools ${on_var} - | gzip > ${on}
        # Link input to output for testing
        ln -s !{sam} !{sample}_in.sam
        '''
}