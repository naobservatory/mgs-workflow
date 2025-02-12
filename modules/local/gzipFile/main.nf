// Gzip a plaintext input file
// Currently only used for testing
process GZIP_FILE {
    label "single"
    label "coreutils"
    input:
        tuple val(sample), path(file)
    output:
        tuple val(sample), path("${file}.gz")
    shell:
        '''
        gzip -c !{file} > !{file}.gz
        '''
}

// Gzip a file without sample annotation
// Currently only used for testing
process GZIP_FILE_BARE {
    label "single"
    label "coreutils"
    input:
        tuple path(file)
    output:
        path("${file}.gz")
    shell:
        '''
        gzip -c !{file} > !{file}.gz
        '''
}
