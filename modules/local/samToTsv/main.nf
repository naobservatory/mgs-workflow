process SAM_TO_TSV {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(sam)
        path genbank_metadata_path
        path viral_db_path
    output:
        tuple val(sample), path("${sample}_bowtie2_sam_processed.tsv.gz"), emit: output
        tuple val(sample), path("${sample}_bowtie2_sam_in.tsv.gz"), emit: input
    shell:
        '''
        out=!{sample}_bowtie2_sam_processed.tsv.gz
        meta=!{genbank_metadata_path}
        db=!{viral_db_path}
        infile=!{sam}
        cmd="process_bowtie2_sam.py -s ${infile} -m ${meta} -v ${db} -o ${out}"
        # Run script
        ${cmd}
        # Link input to output for testing
        ln -s !{sam} !{sample}_bowtie2_sam_in.tsv.gz
        '''
}
