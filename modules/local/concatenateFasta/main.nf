// Concatenate gzipped FASTA files
process CONCATENATE_FASTA_GZIPPED {
    label "single"
    label "base"
    input:
        path(files)
        val(name)
    output:
        path("${name}.fasta.gz")
    shell:
        '''
        cat !{files} > !{name}.fasta.gz
        '''
}

// Concatenate gzipped FASTA files within a directory
process CONCATENATE_FASTA_GZIPPED_DIR {
    label "single"
    label "base"
    input:
        path(dir)
        val(name)
        val(suffix)
    output:
        path("${name}.fasta.gz")
    shell:
        '''
        cat !{dir}/*.!{suffix} > !{name}.fasta.gz
        '''
}
// Concatenate gzipped FASTA files within a directory of subdirectories
process CONCATENATE_FASTA_GZIPPED_DIR_DEEP {
    label "single"
    label "base"
    input:
        path(dir)
        val(name)
        val(suffix)
    output:
        path("${name}.fasta.gz")
    shell:
        '''
        cat !{dir}/*/*.!{suffix} > !{name}.fasta.gz
        '''
}
