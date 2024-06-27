// Concatenate gzipped FASTA files
process CONCATENATE_FASTA_GZIPPED {
    label "single"
    label "base"
    input:
        path(files)
    output:
        path("${params.name}.fasta.gz")
    shell:
        '''
        cat !{files} > !{params.name}.fasta.gz
        '''
}

// Concatenate gzipped FASTA files within a directory
process CONCATENATE_FASTA_GZIPPED_DIR {
    label "single"
    label "base"
    input:
        path(dir)
    output:
        path("${params.name}.fasta.gz")
    shell:
        '''
        cat !{dir}/*.!{params.suffix} > !{params.name}.fasta.gz
        '''
}
