// Process SAM file (add reference taxid, add clean read information, turn into TSV)
process PROCESS_VIRAL_MINIMAP2_SAM {
    label "pysam_biopython"
    label "single"
    input:
        tuple val(sample), path(virus_sam), path(clean_reads)
        path genbank_metadata_path
        path viral_db_path

    output:
        tuple val(sample), path("${sample}_minimap2_sam_processed.tsv.gz"), emit: output
        tuple val(sample), path("input_${virus_sam}"), emit: input
    shell:
        '''
        out=!{sample}_minimap2_sam_processed.tsv.gz
        metadata=!{genbank_metadata_path}
        virus_db=!{viral_db_path}
        virus_sam=!{virus_sam}
        clean_reads=!{clean_reads}
        process_viral_minimap2_sam.py -a ${virus_sam} -r ${clean_reads} -m ${metadata} -v ${virus_db} -o ${out}
        # Link input to output for testing
        ln -s !{virus_sam} input_!{virus_sam}
        '''
}
