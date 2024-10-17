// Annotate virus taxonomic DB with infection status for target host taxa
process ANNOTATE_VIRUS_INFECTION {
    label "single"
    label "biopython"
    input:
        path(virus_db)
        path(host_db)
        path(infection_db)
        path(nodes_db)
    output:
        path("total-virus-db-annotated.tsv.gz"), emit: db
    shell:
        '''
        annotate-viral-hosts.py !{virus_db} !{host_db} !{infection_db} !{nodes_db} total-virus-db-annotated.tsv.gz
        '''
}
