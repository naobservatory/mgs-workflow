// Download genome reference
process DOWNLOAD_GENOME {
    label "BBTools"
    label "single"
    input:
        val(genome_url)
    output:
        path("${params.name}.fasta.gz"), emit: genome
    shell:
        '''
        path=!{params.name}.fasta.gz
        wget !{genome_url} -O ${path}
        '''
}
