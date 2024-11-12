// Download genome reference
process DOWNLOAD_GENOME {
    label "BBTools"
    label "single"
    input:
        val(genome_url)
        val(name)
    output:
        path("${name}.fasta.gz"), emit: genome
    shell:
        '''
        path=!{name}.fasta.gz
        wget !{genome_url} -O ${path}
        '''
}
