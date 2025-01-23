process FILTLONG {
    label "small"
    label "filtlong"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_filtlong.fastq.gz"), emit: reads
    shell:
        // Filter reads based on length (min 100 bp, max 10000 bp) and mean quality (min 99%, i.e, a Phred score of 20)
        '''
        o=!{sample}_filtlong.fastq.gz
        i=!{reads[0]}
        filtlong --min_length 100 --max_length 10000 --min_mean_q 99 --verbose ${i} | gzip > ${o}
        '''
}