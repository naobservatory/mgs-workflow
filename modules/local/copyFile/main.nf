// Copy a file while retaining channel structure
process COPY_FILE {
    label "single"
    label "coreutils"
    input:
        tuple val(sample), path(file)
        val(outname)
    output:
        tuple val(sample), path("${sample}_${outname}")
    shell:
        '''
        cp !{file} !{sample}_!{outname}
        '''
}

// Copy a file without sample annotation
process COPY_FILE_BARE {
    label "single"
    label "coreutils"
    label "testing_only" // Process is currently only used for testing
    input:
        path(file)
        val(outname)
    output:
        path("${outname}")
    shell:
        '''
        cp !{file} !{outname}
        '''
}
