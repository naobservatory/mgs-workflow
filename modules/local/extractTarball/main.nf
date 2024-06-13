// Extract a Gzipped tarball into a directory
process EXTRACT_TARBALL(
    label "BBTools"
    label "single"
    input:
        path(tarball)
    output:
        path("extracted_tarball")
    shell:
        '''
        tar -xzf !{tarball}
        '''
