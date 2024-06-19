// Download human subset of VirusHostDB
// NB: 9606 is the taxid for humans
process DOWNLOAD_HUMAN_VIRUS_DB {
    label "single"
    label "base"
    input:
        val(virus_host_db_url)
    output:
        path("human-virus-db-raw.tsv")
    shell:
        '''
        curl -sS !{virus_host_db_url} | awk -F'\t' '$8=="9606"{print $1"\t"$2}' | sort -n > human-virus-db-raw.tsv
        '''
}
