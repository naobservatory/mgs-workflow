// Filter viral genome data from ncbi-genome-download based on infection status
process FILTER_VIRAL_GENBANK_METADATA {
    label "single"
    label "pandas"
    input:
        path(metadata_db)
        path(virus_db)
        val(host_taxa)
    output:
        path("${params.name}-metadata.tsv.gz"), emit: db
        path("${params.name}-gids.csv"), emit: gid
    shell:
        '''
        filter-viral-genbank-metadata.py !{metadata_db} !{virus_db} "!{host_taxa}" !{params.name}-metadata.tsv.gz !{params.name}-gids.csv
        '''
}
