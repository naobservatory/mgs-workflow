process DOWNLOAD_BLAST_DB {
    label "BLAST2"
    label "max"
    errorStrategy "terminate"
    output:
        path("${params.db}"), emit: db
    shell:
        '''
        db="!{params.db}"
        mkdir ${db}
        cd ${db}
        update_blastdb.pl --decompress ${db}
        '''
}
