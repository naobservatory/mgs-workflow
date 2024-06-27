// Copy a file to a new location with a custom path
process COPY_FILE {
    label "base"
    label "single"
    input:
        path(file_path)
    output:
        path("${params.outpath}"), emit: file
    shell:
        '''
        cp !{file_path} !{params.outpath}
        '''
}
