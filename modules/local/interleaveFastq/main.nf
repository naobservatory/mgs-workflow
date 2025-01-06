// Interleave paired FASTQ files into a single interleaved file

process INTERLEAVE_FASTQ_SEQTK {
    label "seqtk"
    label "single"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_interleaved.fastq.gz")
    shell:
        '''
        seqtk mergepe !{reads[0]} !{reads[1]} | gzip -c > "!{sample}_interleaved.fastq.gz"
        '''
}

process INTERLEAVE_FASTQ_BBTOOLS {
    label "small"
    label "BBTools"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_interleaved.fastq.gz")
    shell:
        '''
        reformat.sh in1=!{reads[0]} in2=!{reads[1]} out="!{sample}_interleaved.fastq.gz"
        '''
}

process INTERLEAVE_FASTQ_CAT {
    label "single"
    label "coreutils_gzip_gawk"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_interleaved.fastq.gz")
    shell:
        '''
        paste <(cat !{reads[0]} | paste - - - - ) <(cat !{reads[1]} | paste - - - - ) | tr "\t" "\n" > "!{sample}_interleaved.fastq.gz"
        '''
}
