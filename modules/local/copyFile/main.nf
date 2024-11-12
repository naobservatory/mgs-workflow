// Copy a file to a new location with a custom path
process COPY_FILE {
    label "base"
    label "single"
    input:
        path(file_path)
        val(outpath)
    output:
        path("${outpath}"), emit: file
    shell:
        '''
        cp !{file_path} !{outpath}
        '''
}
