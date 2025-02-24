// Pull out original clean reads from cleaned fastq

process PULLOUT_FASTQ {
    label "seqtk"
    label "single"
    input:
        tuple val(sample), path(fastq)

    output:
        tuple val(sample), path("${sample}_reads_subset.fastq.gz"), emit: output
        tuple val(sample), path("${sample}_reads_in.fastq.gz"), emit: input
    shell:
        '''
       
        '''
}
