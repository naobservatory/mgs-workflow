// Extract information from a Minimap2 SAM file based on Genbank download metadata
process PROCESS_VIRAL_MINIMAP2_SAM {
    label "pandas"
    label "single"
    input:
        path sam
        path genbank_metadata_path
        path viral_db_path
        val host_taxon
    output:
        path("minimap2_sam_filtered.tsv.gz"), emit: output
        path("minimap2_sam_filtered.sam"), emit: input
    shell:
        '''
        out= minimap2_sam_filtered.tsv.gz
        metadata=!{genbank_metadata_path}
        virus_db=!{viral_db_path}
        host=!{host_taxon}
        process_viral_minimap2_sam.py -s ${sam} -m ${metadata} -v ${virus_db} -ht ${host} -o ${out}
        # Link input to output for testing
        ln -s !{sam} minimap2_sam_filtered.sam
        '''
}
