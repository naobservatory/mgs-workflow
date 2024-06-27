// Extract a Gzipped tarball into a directory
process EXTRACT_TARBALL{
    label "BBTools"
    label "single"
    input:
        path(tarball)
        val(outdir)
        val(makedir)
    output:
        path(outdir)
    shell:
        '''
        if [[ "!{makedir}" == "true" ]]; then
            mkdir !{outdir}
            tar -xzf !{tarball} -C !{outdir}
        else
            tar -xzf !{tarball}
        fi
        '''
}
