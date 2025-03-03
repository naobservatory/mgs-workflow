// Process SAM file(filter out non-host-taxon identifying reads, add reference taxid and name, turn into TSV)
process PROCESS_VIRAL_MINIMAP2_SAM {
    label "pysam_biopython"
    label "single"
    input:
        path hv_sam
        path clean_reads
        path genbank_metadata_path
        path viral_db_path
        val host_taxon

    output:
        path("minimap2_sam_filtered.tsv"), emit: output
        path("input_${hv_sam}"), emit: input
    shell:
        '''

        clean_reads=!{clean_reads}
        out=minimap2_sam_filtered.tsv
        metadata=!{genbank_metadata_path}
        virus_db=!{viral_db_path}
        host=!{host_taxon}
        process_viral_minimap2_sam.py -s !{hv_sam} -r ${clean_reads} -m ${metadata} -v ${virus_db} -ht ${host} -o ${out}
        # Link input to output for testing
        ln -s !{hv_sam} input_!{hv_sam}
        '''
}
