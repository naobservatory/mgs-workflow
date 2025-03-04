// Replace a string in a file while retaining channel structure
// Currently only used for testing
process REPLACE_STRING {
    label "single"
    label "coreutils"
    input:
        tuple val(label), path(file)
        val(outname)
        val(old_string)
        val(new_string)
    output:
        tuple val(label), path("${label}_${outname}")
    shell:
        '''
        sed 's/!{old_string}/!{new_string}/g' !{file} > !{label}_!{outname}
        '''
}

process REPLACE_STRING_IN_COMPRESSED_FILE {
    label "single"
    label "coreutils"
    input:
        tuple val(label), path(file)
        val(outname)
        val(old_string)
        val(new_string)
    output:
        tuple val(label), path("${label}_${outname}")
    shell:
        '''
        zcat !{file} | sed 's/!{old_string}/!{new_string}/g' | gzip > !{label}_!{outname}
        '''
}