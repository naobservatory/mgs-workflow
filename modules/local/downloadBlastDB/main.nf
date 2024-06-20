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
        ln -s $(which curl) /usr/local/bin/curl
        update_blastdb.pl --source aws --num_threads !{task.cpus} --force --decompress ${db}
        '''
}
