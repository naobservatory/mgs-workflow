// Download entire VirusHostDB
process DOWNLOAD_HUMAN_VIRUS_DB {
    label "single"
    label "base"
    input:
        val(virus_host_db_url)
    output:
        path("virus-host-db.tsv")
    shell:
        '''
        curl -sS !{virus_host_db_url} > virus-host-db.tsv
        '''
}
