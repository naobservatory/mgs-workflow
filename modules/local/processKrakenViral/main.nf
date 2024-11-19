// Process Kraken2 output and identify virus- and non-virus-assigned reads
process PROCESS_KRAKEN_VIRAL {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(output)
        path viral_db_path
        val host_taxon
    output:
        tuple val(sample), path("${sample}_kraken_viral_processed.tsv")
    shell:
        '''
        in=!{output}
        out=!{sample}_kraken_viral_processed.tsv
        viral=!{viral_db_path}
        host=!{host_taxon}
        process_kraken_viral.py ${in} ${viral} ${host} ${out}
        '''
}
