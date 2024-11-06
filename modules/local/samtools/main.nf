process BAM_TO_FASTQ {
    label "samtools"

    input:
        path(bam_dir)
        val nanopore_run
    output:
        path '*.fastq.gz'

    shell:
        '''
        # Run samtools
        for bam in !{bam_dir}/*.bam; do
            basename=$(basename "$bam" .bam)
            barcode=$(echo $basename | sed 's/b7f847d7a590c4991a770d9fe21324ef21b88a6c_//')
            samtools fastq "$bam" | gzip -c > "!{nanopore_run}_${barcode}.fastq.gz"        done
        '''
}