// Convert BAM to FASTQ
process BAM_TO_FASTQ {
    label "samtools"

    input:
        path(bam_file)
        val nanopore_run
    output:
        path "${nanopore_run}_${bam_file.baseName}.fastq.gz", emit: fastq

    shell:
        '''
        # Run samtools
        samtools fastq !{bam_file} | gzip -c > !{nanopore_run}_!{bam_file.baseName}.fastq.gz
        '''
}