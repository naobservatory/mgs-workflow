// Process SAM file(filter out non-host-taxon identifying reads, add reference taxid and name, turn into TSV)
process PROCESS_VIRAL_MINIMAP2_SAM {
    label "pysam_biopython"
    label "single"
    input:
        tuple val(sample), path(hv_sam), path(clean_reads)
        path genbank_metadata_path
        path viral_db_path
        val host_taxon

    output:
        tuple val(sample), path("${sample}_minimap2_sam_processed.tsv.gz"), emit: output
        tuple val(sample), path("input_${hv_sam}"), emit: input
    shell:
        '''
        out=!{sample}_minimap2_sam_processed.tsv.gz
        metadata=!{genbank_metadata_path}
        virus_db=!{viral_db_path}
        hv_sam=!{hv_sam}
        clean_reads=!{clean_reads}
        host_taxon=!{host_taxon}
        process_viral_minimap2_sam.py -a ${hv_sam} -r ${clean_reads} -m ${metadata} -v ${virus_db} -t ${host_taxon} -o ${out}
        # Link input to output for testing
        ln -s !{hv_sam} input_!{hv_sam}
        '''
}
