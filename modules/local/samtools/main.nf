// Return reads that did not align to reference as FASTQ
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
