// Combine concordantly and non-concordantly mapped read pairs
process COMBINE_MAPPED_BOWTIE2_READS {
    label "BBTools"
    label "single"
    input:
        tuple val(sample), path(reads_conc), path(reads_mapped_unconc)
    output:
        tuple val(sample), path("${sample}_bowtie2_mapped_all_{1,2}.fastq.gz")
    shell:
        '''
        inc1=!{reads_conc[0]}
        inc2=!{reads_conc[1]}
        inu1=!{reads_mapped_unconc[0]}
        inu2=!{reads_mapped_unconc[1]}
        out1=!{sample}_bowtie2_mapped_all_1.fastq.gz
        out2=!{sample}_bowtie2_mapped_all_2.fastq.gz
        cat $inc1 $inu1 > $out1
        cat $inc2 $inu2 > $out2
        '''
}
