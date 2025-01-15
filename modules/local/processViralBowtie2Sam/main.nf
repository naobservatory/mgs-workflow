// Extract information from a paired Bowtie2 SAM file based on Genbank download metadata

process PROCESS_VIRAL_BOWTIE2_SAM {
    label "pandas"
    label "single"
    input:
        tuple val(sample), path(sam)
        path genbank_metadata_path
        path viral_db_path
    output:
        tuple val(sample), path("${sample}_bowtie2_sam_processed.tsv")
    shell:
        '''
        in=!{sam}
        out=!{sample}_bowtie2_sam_processed.tsv
        meta=!{genbank_metadata_path}
        db=!{viral_db_path}
        process_viral_bowtie2_sam.py ${in} ${meta} ${db} ${out}
        '''
}

// Updated version with testing and gzipped output
process PROCESS_VIRAL_BOWTIE2_SAM_2 {
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
        cmd="process_viral_bowtie2_sam_2.py -m ${meta} -v ${db} -o ${out}"
        # Sort input SAM and pass to script
        zcat !{sam} | sort -t $'\t' -k1,1 | ${cmd}
        # Link input to output for testing
        ln -s !{sam} !{sample}_bowtie2_sam_in.tsv.gz
        '''
}
