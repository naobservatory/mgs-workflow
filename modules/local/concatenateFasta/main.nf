// Concatenate gzipped FASTA files
process CONCATENATE_FASTA_GZIPPED {
    input:
        path(files)
    output:
        path("${params.name}.fasta.gz")
    shell:
        '''
        cat ${files} > !{params.name}.fasta.gz
        '''
}
