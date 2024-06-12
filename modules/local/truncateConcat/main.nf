// Truncate concatenated read files for trial run
process TRUNCATE_CONCAT {
    label "single"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
        val n_reads
    output:
        tuple val(sample), path("${sample}_trunc_{1,2}.fastq.gz")
    shell:
        '''
        echo "Number of output reads: !{n_reads}"
        n_lines=$(expr !{n_reads} \\* 4)
        echo "Number of output lines: ${n_lines}"
        o1=!{sample}_trunc_1.fastq.gz
        o2=!{sample}_trunc_2.fastq.gz
        zcat !{reads[0]} | head -n ${n_lines} | gzip -c > ${o1}
        zcat !{reads[1]} | head -n ${n_lines} | gzip -c > ${o2}
        '''
}
