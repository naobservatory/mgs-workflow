// Copy a file while retaining channel structure (useful for testing)
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
