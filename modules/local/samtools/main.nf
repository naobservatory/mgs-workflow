// Throw away reads that aligned to reference
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
        samtools ${var} ${in} > ${out}
        '''
}

// Return matching and non-matching reads separately
process SAMTOOLS_SEPARATE {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}_samtools_fail.fastq.gz"), emit: fail
        tuple val(sample), path("${sample}_${suffix}_samtools_pass.fastq.gz"), emit: reads
    shell:
        '''
        in=!{sam}
        of=!{sample}_!{suffix}_samtools_fail.fastq.gz
        op=!{sample}_!{suffix}_samtools_pass.fastq.gz
        of_var="fastq -n -F 4"
        op_var="fastq -n -f 4"
        samtools ${of_var} ${in} > ${of}
        samtools ${op_var} ${in} > ${op}
        '''
}