// Process Kraken2 output and identify virus- and non-virus-assigned reads
// Updated version with testing and gzipped output
process PROCESS_KRAKEN_VIRAL {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(output)
        path viral_db_path
        val host_taxon
    output:
        tuple val(sample), path("${sample}_kraken_viral_processed.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_kraken_viral_in.tsv.gz"), emit: input
    shell:
        '''
        in=!{output}
        out=!{sample}_kraken_viral_processed.tsv.gz
        viral=!{viral_db_path}
        host=!{host_taxon}
        process_kraken_viral.py ${in} ${viral} ${host} ${out}
        # Link input for testing
        ln -s !{output} !{sample}_kraken_viral_in.tsv.gz
        '''
}
