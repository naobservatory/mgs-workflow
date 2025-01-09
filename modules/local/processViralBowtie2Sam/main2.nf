// Extract information from a paired Bowtie2 SAM file based on Genbank download metadata

process PROCESS_VIRAL_BOWTIE2_SAM {
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
        in=!{sam}
        out=!{sample}_bowtie2_sam_processed.tsv.gz
        meta=!{genbank_metadata_path}
        db=!{viral_db_path}
        process_viral_bowtie2_sam.py ${in} ${meta} ${db} ${out}
        ln -s !{sam} !{sample}_bowtie2_sam_in.tsv.gz
        '''
}
