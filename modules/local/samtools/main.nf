process BAM_TO_FASTQ {
    label "samtools"

    input:
        path(bam)
        val nanopore_run
    output:
        path '*.fastq.gz'

    shell:
        '''
        base_name=$(basename !{bam} .bam)
        # Remove "calls_" prefix from base_name in case of no demultiplexing
        base_name=${base_name#calls_}
        samtools fastq !{bam} | gzip -c > "${base_name}.fastq.gz"
        '''
}