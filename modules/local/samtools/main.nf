// Return reads unaligned reads as FASTQ
process SAMTOOLS_FILTER {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}.fastq.gz"), emit: reads
    shell:
        '''
        in=!{sam}
        out=!{sample}_!{suffix}.fastq.gz
        var="fastq -n -f 4"
        samtools ${var} ${in} | gzip > ${out}
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
    shell:
        '''
        in=!{sam}
        of=!{sample}_!{suffix}_samtools_fail.fastq.gz
        op=!{sample}_!{suffix}_samtools_pass.fastq.gz
        of_var="fastq -n -F 4"
        op_var="fastq -n -f 4"
        samtools ${of_var} ${in} | gzip > ${of}
        samtools ${op_var} ${in} | gzip > ${op}
        '''
}