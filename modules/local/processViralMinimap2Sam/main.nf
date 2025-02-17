// Extract information from a Minimap2 SAM file based on Genbank download metadata
process PROCESS_VIRAL_MINIMAP2_SAM {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(sam)
        path genbank_metadata_path
        path viral_db_path
        path host_taxon
    output:
        tuple val(sample), path("${sample}_minimap2_sam_processed.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_minimap2_in.sam"), emit: input
    shell:
        '''
        out=!{sample}_minimap2_sam_processed.tsv.gz
        metadata=!{genbank_metadata_path}
        virus_db=!{viral_db_path}
        host=!{host_taxon}
        process_viral_minimap2_sam.py -s ${sam} -m ${metadata} -v ${virus_db} -ht ${host} -o ${out}
        # Link input to output for testing
        ln -s !{sam} !{sample}_minimap2_in.sam
        '''
}
