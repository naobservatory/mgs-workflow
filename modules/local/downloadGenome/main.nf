// Download genome reference
process DOWNLOAD_GENOME {
    label "BBTools"
    label "single"
    input:
        tuple val(genome_url), val(name)
    output:
        path("${name}.fasta.gz"), emit: genome
    shell:
        '''
        path=!{name}.fasta
        wget "!{genome_url}" -O ${path}
        gzip ${path}
        '''
}
