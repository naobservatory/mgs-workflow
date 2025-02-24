// Process SAM file(filter out non-host-taxon identifying reads, add reference taxid and name, turn into TSV)
process PROCESS_VIRAL_MINIMAP2_SAM {
    label "pysam"
    label "single"
    input:
        path sam
        path genbank_metadata_path
        path viral_db_path
        val host_taxon
    output:
        path("minimap2_sam_filtered.tsv"), emit: output
        path("minimap2_sam_filtered.sam"), emit: input
    shell:
        '''
        in=!{sam}
        out=minimap2_sam_filtered.tsv
        metadata=!{genbank_metadata_path}
        virus_db=!{viral_db_path}
        host=!{host_taxon}
        process_viral_minimap2_sam.py -s ${in} -m ${metadata} -v ${virus_db} -ht ${host} -o ${out}
        # Link input to output for testing
        ln -s !{sam} minimap2_sam_filtered.sam
        '''
}
