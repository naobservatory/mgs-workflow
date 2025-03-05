// Concatenate gzipped FASTA files
process CONCATENATE_FASTQ_GZIPPED {
    label "single"
    label "coreutils"
    input:
        path(files)
        val(name)
    output:
        path("${name}.fastq.gz")
    shell:
        '''
        cat !{files} > !{name}.fastq.gz
        '''
}