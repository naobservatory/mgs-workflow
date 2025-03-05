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
        sample=!{sample}
        metadata=!{genbank_metadata_path}
        virus_db=!{viral_db_path}
        hv_sam=!{hv_sam}
        clean_reads=!{clean_reads}
        host=!{host_taxon}
        process_viral_minimap2_sam.py -sa ${hv_sam} -r ${clean_reads} -sl ${sample} -m ${metadata} -v ${virus_db} -ht ${host} -o ${out}
        # Link input to output for testing
        ln -s !{hv_sam} input_!{hv_sam}
        '''
}
