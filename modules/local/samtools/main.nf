process SAMTOOLS_FILTER_CONTAM {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}.fastq.gz"), emit: reads
    shell:
        '''
        in=!{sam}
        out=!{sample}_${suffix}.fastq.gz
        var="fastq -n -f 4"
        samtools ${var} ${in} -> ${out}
        '''
}

process SAMTOOLS_KEEP_CONTAM {
    label "samtools"
    input:
        tuple val(sample), path(sam)
        val(suffix)
    output:
        tuple val(sample), path("${sample}_${suffix}.fastq.gz"), emit: reads
    shell:
        '''
        in=!{sam}
        out=!{sample}_${suffix}.fastq.gz
        var="fastq -n -F 4"
        samtools ${var} ${in} -> ${out}
        '''

}

