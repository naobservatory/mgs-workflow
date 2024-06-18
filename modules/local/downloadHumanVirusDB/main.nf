// Download human subset of VirusHostDB
process DOWNLOAD_HUMAN_VIRUS_DB {
    input:
        path(virus_host_db_url)
    output:
        path("human-virus-db-raw.tsv")
    shell:
        '''
        curl -sS !{virus_host_db_url} \
        | awk -F'\t' '$8=="9606"{print $1"\t"$2}' \ # 9606 is the taxid for humans
        | sort -n \
        > human-virus-db-raw.tsv
        '''
}
