// Join and concatenate partially-merged FASTQ files into a single read file
process JOIN_FASTQ {
    label "biopython"
    label "single"
    input:
        tuple val(sample), path(reads)
    output:
        tuple val(sample), path("${sample}_joined.fastq.gz")
    shell:
        '''
        # Prepare to join unmerged read pairs
        om=!{reads[0]}
        ou1=!{reads[1]}
        ou2=!{reads[2]}
        oj=!{sample}_bbmerge_unmerged_joined.fastq.gz
        # Join unmerged read pairs
        join_fastq.py ${ou1} ${ou2} ${oj}
        # Concatenate single output file
        oo=!{sample}_joined.fastq.gz
        cat ${om} ${oj} > ${oo}
        '''
}
